/**
 * @file    codegen_cc.cc
 * @date    201810xx
 * @author  Wouter Vlothuizen (wouter.vlothuizen@tno.nl)
 * @brief   code generator backend for the Central Controller
 * @note    here we don't check whether the sequence of calling code generator
 *          functions is correct
 */

#include "codegen_cc.h"

#include <version.h>
#include <options.h>
#include <iosfwd>

// constants


namespace ql {


#if OPT_FEEDBACK
//typedef std::vector<int> tBitVars;
#endif


/************************************************************************\
| Generic
\************************************************************************/

void codegen_cc::init(const quantum_platform &platform)
{
    // NB: a new eqasm_backend_cc is instantiated per call to compile, and
    // as a result also a codegen_cc, so we don't need to cleanup
    this->platform = &platform;
    settings.loadBackendSettings(platform);

    // optionally preload codewordTable
    std::string map_input_file = options::get("backend_cc_map_input_file");
    if(!map_input_file.empty()) {
        QL_DOUT("loading map_input_file='" << map_input_file << "'");
        Json map = load_json(map_input_file);
        codewordTable = map["codeword_table"];      // FIXME: use json_get
        mapPreloaded = true;
    }

#if OPT_FEEDBACK   // FIXME: WIP on feedback: allocate SM bits
    // iterate over instruments
    for(size_t instrIdx=0; instrIdx<settings.getInstrumentsSize(); instrIdx++) {
        const settings_cc::tInstrumentControl ic = settings.getInstrumentControl(instrIdx);
        if(QL_JSON_EXISTS(ic.controlMode, "result_bits")) {  // this instrument mode produces results (i.e. it is a measurement device)
            // FIXME: maintain mapping instrument -> SM
            QL_IOUT("instrument '" << ic.ii.instrumentName << "' (index " << instrIdx << ") is used for feedback");
        }
    }
#endif
}

std::string codegen_cc::getProgram()
{
#if OPT_FEEDBACK
	return codeSection.str() + dp.getDatapathSection();
#else
    return codeSection.str();
#endif
}

std::string codegen_cc::getMap()
{
    Json map;

    map["note"] = "generated by OpenQL CC backend version " CC_BACKEND_VERSION_STRING;
    map["codeword_table"] = codewordTable;
    return QL_SS2S(std::setw(4) << map << std::endl);
}


/************************************************************************\
| 'Program' level functions
\************************************************************************/

void codegen_cc::programStart(const std::string &progName)
{
    // emit program header
    codeSection << std::left;    // assumed by emit()
    codeSection << "# Program: '" << progName << "'" << std::endl;   // NB: put on top so it shows up in internal CC logging
    codeSection << "# CC_BACKEND_VERSION " << CC_BACKEND_VERSION_STRING << std::endl;
    codeSection << "# OPENQL_VERSION " << OPENQL_VERSION_STRING << std::endl;
    codeSection << "# Note:    generated by OpenQL Central Controller backend" << std::endl;
    codeSection << "#" << std::endl;
    emitProgramStart();

	dp.programStart();

    vcd.programStart(platform->qubit_number, platform->cycle_time, MAX_GROUPS, settings);
}


void codegen_cc::programFinish(const std::string &progName)
{
#if OPT_RUN_ONCE   // program runs once only
    emit("", "stop");
#else   // CC-light emulation: loop indefinitely
    emit("",      // no CCIO selector
         "jmp",
         "@mainLoop",
         "# loop indefinitely");
#endif

    emit(".END");   // end .CODE section

	dp.programFinish();

    vcd.programFinish(progName);
}

/************************************************************************\
| 'Kernel' level functions
\************************************************************************/

void codegen_cc::kernelStart()
{
    zero(lastEndCycle);       // FIXME: actually, bundle.startCycle starts counting at 1
}

void codegen_cc::kernelFinish(const std::string &kernelName, size_t durationInCycles)
{
    vcd.kernelFinish(kernelName, durationInCycles);
}

/************************************************************************\
| 'Bundle' level functions
\************************************************************************/

/*
    Our strategy is to first process all customGate's in a bundle, storing the
    relevant information in bundleInfo. Then, when all work for a bundle has
    been collected, we generate code in bundleFinish

    - bundleStart():
    clear bundleInfo, which maintains the work that needs to be performed for bundle

    - customGate():
    collect gate information in bundleInfo

    - bundleFinish():
    generate code for bundle from information collected in bundleInfo (which
    may be empty if no custom gates are present in bundle)
*/

// bundleStart: see 'strategy' above
void codegen_cc::bundleStart(const std::string &cmnt)
{
	bundleHasFeedback = false;

    // create 'matrix' of BundleInfo with proper vector size per instrument
	bundleInfo.clear();
    BundleInfo empty;
    for(size_t instrIdx=0; instrIdx<settings.getInstrumentsSize(); instrIdx++) {
        const settings_cc::tInstrumentControl ic = settings.getInstrumentControl(instrIdx);
        bundleInfo.emplace_back(
        	ic.controlModeGroupCnt,   	// one BundleInfo per group in the control mode selected for instrument
			empty  						// empty BundleInfo
		);
    }

	// generate source code comments
    comment(cmnt);
    dp.comment(cmnt);		// FIXME: comment is not fully appropriate, but at least allows matching with .CODE section
}




// Static helper function for bundleFinish()
typedef struct {
	uint32_t groupDigOut;   // codeword/mask fragment for this group
	std::string comment;    // comment for instruction stream
} tCalcGroupDigOut;

static tCalcGroupDigOut calcGroupDigOut(size_t instrIdx, size_t group, size_t nrGroups, const settings_cc::tInstrumentControl &ic, int staticCodewordOverride)
{
	tCalcGroupDigOut ret{0, ""};

	// determine control mode group FIXME: more explanation
	int controlModeGroup = -1;
	if(ic.controlModeGroupCnt == 0) {
		QL_JSON_FATAL("'control_bits' not defined or empty in 'control_modes/" << ic.refControlMode <<"'");
#if OPT_VECTOR_MODE
		} else if(ic.controlModeGroupCnt == 1) {                    // vector mode: group addresses channel within vector
	controlModeGroup = 0;
#endif
	} else if(group < ic.controlModeGroupCnt) {                 // normal mode: group selects control group
		controlModeGroup = group;
	} else {
		// NB: this actually an error in program logic
		QL_JSON_FATAL(
			"instrument '" << ic.ii.instrumentName
			<< "' uses " << nrGroups
			<< " groups, but control mode '" << ic.refControlMode
			<< "' only defines " << ic.controlModeGroupCnt
			<< " groups in 'control_bits'"
		);
	}

	// get number of control bits for group
	const Json &groupControlBits = ic.controlMode["control_bits"][controlModeGroup];    // NB: tests above guarantee existence
	QL_DOUT(
		"instrumentName=" << ic.ii.instrumentName
		<< ", slot=" << ic.ii.slot
		<< ", control mode group=" << controlModeGroup
		<< ", group control bits: " << groupControlBits
	);
	size_t nrGroupControlBits = groupControlBits.size();


	// calculate digital output for group
	if(nrGroupControlBits == 1) {       // single bit, implying this is a mask (not code word)
		ret.groupDigOut |= 1<<(int)groupControlBits[0];     // NB: we assume the mask is active high, which is correct for VSM and UHF-QC
		// FIXME: check controlModeGroup vs group
	} else if(nrGroupControlBits > 1) {                 // > 1 bit, implying code word
#if OPT_VECTOR_MODE
		//  allow single code word for vector of groups. FIXME: requires looking at all sd.signal before assigning code word
	if(group != controlModeGroup) {
		// FIXME: unfinished work on vector mode
	}
#endif

		// find or assign code word
		uint32_t codeword = 0;
		bool codewordOverriden = false;
#if OPT_SUPPORT_STATIC_CODEWORDS
		codeword = staticCodewordOverride;
		codewordOverriden = true;
#else
		codeword = assignCodeword(ic.ii.instrumentName, instrIdx, group);
#endif

		// convert codeword to digOut
		for(size_t idx=0; idx<nrGroupControlBits; idx++) {
			int codeWordBit = nrGroupControlBits-1-idx;    // NB: groupControlBits defines MSB..LSB
			if(codeword & (1<<codeWordBit)) ret.groupDigOut |= 1<<(int)groupControlBits[idx];
		}

		ret.comment = QL_SS2S(
			"  # slot=" << ic.ii.slot
			<< ", instrument='" << ic.ii.instrumentName << "'"
			<< ", group=" << group
			<< ": codeword=" << codeword
			<< std::string(codewordOverriden ? " (static override)" : "")
			<< ": groupDigOut=0x" << std::hex << std::setfill('0') << std::setw(8) << ret.groupDigOut
		);
	} else {    // nrGroupControlBits < 1
		QL_JSON_FATAL(
			"key 'control_bits' empty for group " << controlModeGroup
			<< " on instrument '" << ic.ii.instrumentName << "'"
		);
	}

	// add trigger to digOut
	size_t nrTriggerBits = ic.controlMode["trigger_bits"].size();
	if(nrTriggerBits == 0) {                                    // no trigger
		// do nothing
	} else if(nrTriggerBits == 1) {                             // single trigger for all groups (NB: will possibly assigned multiple times)
		ret.groupDigOut |= 1 << (int)ic.controlMode["trigger_bits"][0];
#if 1	// FIXME: hotfix for QWG, implement properly
	} else if(nrTriggerBits == 2) {
        ret.groupDigOut |= 1 << (int)ic.controlMode["trigger_bits"][0];
        ret.groupDigOut |= 1 << (int)ic.controlMode["trigger_bits"][1];
#endif
#if 1   // FIXME: trigger per group
	} else if(nrTriggerBits == nrGroups) {                      // trigger per group
		ret.groupDigOut |= 1 << (int)ic.controlMode["trigger_bits"][group];
#endif
	} else {
		QL_JSON_FATAL(
			"instrument '" << ic.ii.instrumentName
			<< "' uses " << nrGroups
			<< " groups, but control mode '" << ic.refControlMode
			<< "' defines " << nrTriggerBits
			<< " trigger bits in 'trigger_bits' (must be 1 or #groups)"
		);
	}

	return ret;
}


// bundleFinish: see 'strategy' above
// FIXME: split into smaller parts
void codegen_cc::bundleFinish(size_t startCycle, size_t durationInCycles, bool isLastBundle)
{
#if OPT_PRAGMA
	bool bundleHasPragma = false;
#endif

    // iterate over instruments
    for(size_t instrIdx=0; instrIdx<settings.getInstrumentsSize(); instrIdx++) {
    	// get control info from instrument settings
        const settings_cc::tInstrumentControl ic = settings.getInstrumentControl(instrIdx);
        if(ic.ii.slot >= MAX_SLOTS) {
            QL_JSON_FATAL(
            	"illegal slot " << ic.ii.slot
            	<< " on instrument '" << ic.ii.instrumentName
			);
        }

		/************************************************************************\
		| collect code generation info from all groups within one instrument
		\************************************************************************/

        // FIXME: the term 'group' is used in a diffused way: 1) index of signal vectors, 2) controlModeGroup
        bool instrHasOutput = false;
        uint32_t digOut = 0;                                                // the digital output value sent over the instrument interface
        unsigned int instrMaxDurationInCycles = 0;                          // maximum duration over groups that are used, one instrument
#if OPT_FEEDBACK
		tFeedbackMap feedbackMap;
		tCondGateMap condGateMap;
#endif
#if OPT_PRAGMA
		const Json *pragma = nullptr;
		int pragmaSmBit = 0;
#endif

		// now collect code generation info from all groups of instrument
        size_t nrGroups = bundleInfo[instrIdx].size();
        for(size_t group=0; group<nrGroups; group++) {
            BundleInfo *bi = &bundleInfo[instrIdx][group];                 	// shorthand

            // handle output
            if(!bi->signalValue.empty()) {                                  // signal defined, i.e.: we need to output something
                // compute maximum duration over all groups
                if(bi->durationInCycles > instrMaxDurationInCycles) instrMaxDurationInCycles = bi->durationInCycles;

                tCalcGroupDigOut gdo = calcGroupDigOut(instrIdx, group, nrGroups, ic, bi->staticCodewordOverride);
                digOut |= gdo.groupDigOut;
                comment(gdo.comment);
#if OPT_FEEDBACK
				// conditional gates
				// store condition and groupDigOut in condMap, if all groups are unconditional we use old scheme, otherwise
				// datapath is configured to generate proper digital output
                if(bi->condition==cond_always || ic.ii.forceCondGatesOn) {
                	// nothing to do, just use digOut
                } else {	// other conditions, including cond_never
					// remind mapping for setting PL
					condGateMap.emplace(group, tCondGateInfo{bi->condition, bi->cond_operands, gdo.groupDigOut});
                }
#endif

                vcd.bundleFinishGroup(startCycle, bi->durationInCycles, gdo.groupDigOut, bi->signalValue, instrIdx, group);

                instrHasOutput = true;
            } // if(signal defined)


#if OPT_PRAGMA
			// handle pragma
			if(bi->pragma) {
				// FIXME: enforce single pragma per bundle (currently by design)
				// FIXME: enforce no other work
            	bundleHasPragma = true;
            	pragma = bi->pragma;

				// FIXME: use breg_operands if present? How about qubit (operand) then?
				int breg_operand = bi->operands[0];                	// implicit classic bit for qubit. FIXME: perform checks
				// get SM bit for classic operand (allocated during readout)
				pragmaSmBit = dp.getSmBit(breg_operand, instrIdx);
            }
#endif


#if OPT_FEEDBACK
            // handle readout (i.e. when necessary, create feedbackMap entry
            // NB: we allow for instruments that perform the input side of readout only, without signal generation by the
            // same instrument, which might be needed in the future
            // FIXME: also generate VCD

			if(bi->isMeasFeedback) {
				int resultBit = settings_cc::getResultBit(ic, group);

#if 0	// FIXME: partly redundant
				// get our qubit
				const Json qubits = json_get<const Json>(*ic.ii.instrument, "qubits", ic.ii.instrumentName);   // NB: json_get<const Json&> unavailable
				size_t qubitGroupCnt = qubits.size();                                  // NB: JSON key qubits is a 'matrix' of [groups*qubits]
				if (group >= qubitGroupCnt) {	// FIXME: also tested in settings_cc::findSignalInfoForQubit
					QL_FATAL("group " << group << " not defined in '" << ic.ii.instrumentName << "/qubits'");
				}
				const Json qubitsOfGroup = qubits[group];
				if (qubitsOfGroup.size() != 1) {	// FIXME: not tested elsewhere
					QL_FATAL("group " << group << " of '" << ic.ii.instrumentName << "/qubits' should define 1 qubit, not " << qubitsOfGroup.size());
				}
				int qubit = qubitsOfGroup[0];
				if (bi->readoutQubit != qubit) {          	// this instrument group handles requested qubit. FIXM: inherently true
					QL_FATAL("inconsistency FIXME");
				};
#endif
				// get classic operand
				if (!bi->breg_operands.empty()) {	// FIXME: breg_operands should be honoured
					// FIXME: Note that our gate decomposition "measure_fb %0": ["measure %0", "_wait_uhfqa %0", "_dist_dsm %0", "_wait_dsm %0"] is problematic in that sense
					QL_WOUT("ignoring explicit assignment to bit " << bi->breg_operands[0] << " for measurement of qubit " << bi->operands[0]);
				}
				int breg_operand = bi->operands[0];                	// implicit classic bit for qubit

				// allocate SM bit for classic operand
				int smBit = dp.allocateSmBit(breg_operand, instrIdx);

				// remind mapping of bit -> smBit for setting MUX
				feedbackMap.emplace(group, tFeedbackInfo{smBit, resultBit, bi});
			}
#endif
        } // for(group)


		/************************************************************************\
		| turn code generation info collected above into actual code
		\************************************************************************/

		if(isLastBundle && instrIdx==0) {
			comment(QL_SS2S(" # last bundle of kernel, will pad outputs to match durations"));
		}

        // generate code for instrument output
        if(instrHasOutput) {
			emitOutput(condGateMap, digOut, instrMaxDurationInCycles, instrIdx, startCycle, ic.ii.slot, ic.ii.instrumentName);
		} else {    // !instrHasOutput
			// nothing to do, we delay emitting till a slot is used or kernel finishes (i.e. isLastBundle just below)
		}

#if OPT_PRAGMA
		if(pragma) {	// NB: note that this will only work because we set the pragma for all instruments, and thus already encounter this for the first instrument
			emitPragma(pragma, pragmaSmBit, instrIdx, startCycle, ic.ii.slot, ic.ii.instrumentName);
        }
#endif

#if OPT_FEEDBACK
		if(bundleHasFeedback) {
			emitFeedback(feedbackMap, instrIdx, startCycle, ic.ii.slot, ic.ii.instrumentName);
		}
#endif

		// for last bundle, pad end of bundle to align durations
        if(isLastBundle) {
            padToCycle(instrIdx, startCycle+durationInCycles, ic.ii.slot, ic.ii.instrumentName);		// FIXME: use instrMaxDurationInCycles and/or check consistency
        }

        vcd.bundleFinish(startCycle, digOut, instrMaxDurationInCycles, instrIdx);	// FIXME: conditional gates, etc
    } // for(instrIdx)

    comment("");    // blank line to separate bundles
}

/************************************************************************\
| Quantum instructions
\************************************************************************/

// helper
static std::string qasm(const std::string &iname, const Vec<UInt> &operands, const Vec<UInt> &breg_operands)
{
	// FIXME: hack
	custom_gate g(iname);
	g.operands = operands;
	g.breg_operands = breg_operands;
	return g.qasm();
}

// customGate: single/two/N qubit gate, including readout, see 'strategy' above
// translates 'gate' representation to 'waveform' representation (BundleInfo) and maps qubits to instruments & group.
// Does not deal with the control mode and digital interface of the instrument.
void codegen_cc::customGate(
	const std::string &iname,
	const Vec<UInt> &operands,
	const Vec<UInt> &creg_operands,
	const Vec<UInt> &breg_operands,
	cond_type_t condition,
	const Vec<UInt> &cond_operands,
	double angle,
	size_t startCycle, size_t durationInCycles
)
{
#if 0   // FIXME: test for angle parameter
    if(angle != 0.0) {
        QL_DOUT("iname=" << iname << ", angle=" << angle);
    }
#endif

    vcd.customGate(iname, operands, startCycle, durationInCycles);

    bool isReadout = settings.isReadout(iname);    	//  determine whether this is a readout instruction

    // generate comment
    if(isReadout) {
		comment(std::string(" # READOUT: '") + qasm(iname, operands, breg_operands) + "'");
    } else { // handle all other instruction types than "readout"
        // generate comment. NB: we don't have a particular limit for the number of operands
        comment(std::string(" # gate '") + qasm(iname, operands, breg_operands) + "'");
    }


    // find instruction (gate definition)
    const Json &instruction = platform->find_instruction(iname);
    // find signal vector definition for instruction
    settings_cc::tSignalDef sd = settings.findSignalDefinition(instruction, iname);

    // scatter signals defined for instruction (e.g. several operands and/or types) to instruments & groups
    for(size_t s=0; s<sd.signal.size(); s++) {
        tCalcSignalValue csv = calcSignalValue(sd, s, operands, iname);

        // store signal value, checking for conflicts
        BundleInfo *bi = &bundleInfo[csv.si.instrIdx][csv.si.group];       	// shorthand
        if(!csv.signalValueString.empty()) {								// empty implies no signal
			if(bi->signalValue.empty()) {                                   // signal not yet used
				bi->signalValue = csv.signalValueString;
#if OPT_SUPPORT_STATIC_CODEWORDS
				// FIXME: this does not only provide support, but findStaticCodewordOverride() currently actually requires static codewords
				bi->staticCodewordOverride = ql::settings_cc::findStaticCodewordOverride(instruction, csv.operandIdx, iname); // NB: function return -1 means 'no override'
#endif
			} else if(bi->signalValue == csv.signalValueString) {           // signal unchanged
				// do nothing
			} else {
				showCodeSoFar();
				QL_FATAL(
					"Signal conflict on instrument='" << csv.si.ic.ii.instrumentName
					<< "', group=" << csv.si.group
					<< ", between '" << bi->signalValue
					<< "' and '" << csv.signalValueString << "'"
				);  // FIXME: add offending instruction
			}
		}

        // store signal duration
        bi->durationInCycles = durationInCycles;

#if OPT_FEEDBACK
        // FIXME: assumes that group configuration for readout input matches that of output
		// store operands used for readout, actual work is postponed to bundleFinish()
        if(isReadout) {
			/*
			 * kernel->gate allows 3 types of measurement:
			 * 		- no explicit result. Historically this implies either:
			 * 			- no result, measurement results are often read offline from the readout device (mostly the raw values
			 * 			instead of the binary result), without the control device ever taking notice of the value
			 * 			- implicit bit result for qubit, e.g. for the CC-light using conditional gates
			 * 		- creg result (old. FIXME: what are intended semantics?)
			 * 			note that Creg's are managed through a class, whereas bregs are just numbers
			 * 		- breg result (new)
			 */

			// operand checks
			if(operands.size() != 1) {
				QL_FATAL(
					"Readout instruction '" << qasm(iname, operands, breg_operands)
											<< "' requires exactly 1 quantum operand, not " << operands.size()
				);
			}
			if(!creg_operands.empty()) {
				QL_FATAL("Using Creg as measurement target is deprecated, use new bit registers");
			}
			if(breg_operands.size() > 1) {
				QL_FATAL(
					"Readout instruction '" << qasm(iname, operands, breg_operands)
											<< "' requires 0 or 1 bit operands, not " << breg_operands.size()
				);
			}

			// store operands
			if(settings.getReadoutMode(iname)=="feedback") {
				bundleHasFeedback = true;
				bi->isMeasFeedback = true;
				bi->operands = operands;
//            	bi->creg_operands = creg_operands;	// NB: will be empty because of checks performed earlier
				bi->breg_operands = breg_operands;
			}
        }

        // store 'expression' for conditional gates
        bi->condition = condition;
        bi->cond_operands = cond_operands;
#endif

        QL_DOUT("customGate(): iname='" << iname <<
             "', duration=" << durationInCycles <<
             " [cycles], instrIdx=" << csv.si.instrIdx <<
             ", group=" << csv.si.group);

        // NB: code is generated in bundleFinish()
    }   // for(signal)

#if OPT_PRAGMA
	const Json *pragma = settings.getPragma(iname);
	if(pragma) {
		for(std::vector<BundleInfo> &vbi : bundleInfo) {
			// FIXME: for now we just store the JSON of the pragma statement in bundleInfo[*][0]
			if(vbi[0].pragma) {
				QL_FATAL("Bundle contains more than one gate with 'pragma' key");	// FIXME: provide context
			}
			vbi[0].pragma = pragma;

			// store operands
			vbi[0].operands = operands;
//			vbi[0].creg_operands = creg_operands;	// NB: will be empty because of checks performed earlier
			vbi[0].breg_operands = breg_operands;
		}
	}
#endif
}

void codegen_cc::nopGate()
{
    comment("# NOP gate");
    QL_FATAL("FIXME: NOP gate not implemented");
}

/************************************************************************\
| Classical operations on kernels
\************************************************************************/

void codegen_cc::ifStart(size_t op0, const std::string &opName, size_t op1)
{
    comment(QL_SS2S("# IF_START(R" << op0 << " " << opName << " R" << op1 << ")"));
    QL_FATAL("FIXME: not implemented");
}

void codegen_cc::elseStart(size_t op0, const std::string &opName, size_t op1)
{
    comment(QL_SS2S("# ELSE_START(R" << op0 << " " << opName << " R" << op1 << ")"));
    QL_FATAL("FIXME: not implemented");
}

void codegen_cc::forStart(const std::string &label, int iterations)
{
    comment(QL_SS2S("# FOR_START(" << iterations << ")"));
    // FIXME: reserve register
    emit("", "move", QL_SS2S(iterations << ",R62"), "# R62 is the 'for loop counter'");        // FIXME: fixed reg, no nested loops
    emit((label+":"), "", "", "# ");        // just a label
#if OPT_PRAGMA
    pragmaForLabel = label;		// remind label for pragma/break FIXME: implement properly later on
#endif
}

void codegen_cc::forEnd(const std::string &label)
{
    comment("# FOR_END");
    // FIXME: free register
    emit("", "loop", QL_SS2S("R62,@" << label), "# R62 is the 'for loop counter'");        // FIXME: fixed reg, no nested loops
#if OPT_PRAGMA
    emit((label+"_end:"), "", "", "# ");                          // NB: just a label
#endif
}

void codegen_cc::doWhileStart(const std::string &label)
{
    comment("# DO_WHILE_START");
    emit((label+":"), "", "", "# ");                              // NB: just a label
}

void codegen_cc::doWhileEnd(const std::string &label, size_t op0, const std::string &opName, size_t op1)
{
    comment(QL_SS2S("# DO_WHILE_END(R" << op0 << " " << opName << " R" << op1 << ")"));
    emit("", "jmp", QL_SS2S("@" << label), "# FIXME: we don't support conditions, just an endless loop'");        // FIXME: just endless loop
}

void codegen_cc::comment(const std::string &c)
{
    if(verboseCode) emit(c);
}

/************************************************************************\
|
| private functions
|
\************************************************************************/

/************************************************************************\
| Some helpers to ease nice assembly formatting
\************************************************************************/

// FIXME: assure space between fields!
// FIXME: make comment output depend on verboseCode

void codegen_cc::emit(const std::string &labelOrComment, const std::string &instr)
{
    if(labelOrComment.length()==0) {  					// no label
        codeSection << "        " << instr << std::endl;
    } else if(labelOrComment.length()<8) {              // label fits before instr
        codeSection << std::setw(8) << labelOrComment << instr << std::endl;
    } else if(instr.length()==0) {                      // no instr
        codeSection << labelOrComment << std::endl;
    } else {
        codeSection << labelOrComment << std::endl << "        " << instr << std::endl;
    }
}


// @param   labelOrSel      label must include trailing ":"
// @param   comment     	must include leading "#"
void codegen_cc::emit(const std::string &labelOrSel, const std::string &instr, const std::string &ops, const std::string &comment)
{
    codeSection << std::setw(16) << labelOrSel << std::setw(16) << instr << std::setw(24) << ops << comment << std::endl;
}

void codegen_cc::emit(int sel, const std::string &instr, const std::string &ops, const std::string &comment)
{
	emit(QL_SS2S("[" << sel << "]"), instr, ops, comment);
}

/************************************************************************\
| helpers
\************************************************************************/

void codegen_cc::emitProgramStart()
{
#if OPT_FEEDBACK
    emit(".CODE");   // start .CODE section
#endif

	// NB: new seq_bar semantics (firmware from 20191219 onwards)
    comment("# synchronous start and latency compensation");
    emit("",                "seq_bar",  "",                 "# synchronization, delay set externally through SET_SEQ_BAR_CNT");

    emit("mainLoop:",       "",         "",                 "# ");

#if OPT_FEEDBACK
    // initialize state
    emit("",                "seq_state","0",                "# clear Programmable Logic state");
#endif
}


// generate code to input measurement results and distribute them via DSM
void codegen_cc::emitFeedback(const tFeedbackMap &feedbackMap, size_t instrIdx, size_t startCycle, int slot, const std::string &instrumentName)
{
	if(startCycle > lastEndCycle[instrIdx]) {	// i.e. if(!instrHasOutput)
		padToCycle(instrIdx, startCycle, slot, instrumentName);
	}

	// code generation for participating and non-participating instruments (NB: must take equal number of sequencer cycles)
	if(!feedbackMap.empty()) {	// this instrument performs readout for feedback now
		int smAddr = 0;		// FIXME:
		int mux = dp.getOrAssignMux(instrIdx, feedbackMap);

		dp.emitMux(mux, smAddr, feedbackMap, instrIdx, slot);

		// emit code for slot input
		int sizeTag = datapath_cc::getSizeTag(feedbackMap.size());		// compute DSM transfer size tag (for 'seq_in_sm' instruction)
		emit(
			slot,
			"seq_in_sm",
			QL_SS2S("S" << smAddr << ","  << mux << "," << sizeTag),
			QL_SS2S("# cycle " << lastEndCycle[instrIdx] << "-" << lastEndCycle[instrIdx]+1 << ": feedback on '" << instrumentName+"'")
		);
		lastEndCycle[instrIdx]++;
	} else {	// this instrument does not perform readout for feedback now
		// FIXME:
		int smAddr = 0;
		int smTotalSize = 6;	// FIXME: calculate, requires overview over all measurements of bundle, or take a safe max

		// emit code for non-participating instrument
		emit(
			slot,
			"seq_inv_sm",
			QL_SS2S("S" << smAddr << ","  << smTotalSize),
			QL_SS2S("# cycle " << lastEndCycle[instrIdx] << "-" << lastEndCycle[instrIdx]+1 << ": invalidate SM on '" << instrumentName+"'")
		);
		lastEndCycle[instrIdx]++;
	}
}


void codegen_cc::emitOutput(const tCondGateMap &condGateMap, int32_t digOut, unsigned int instrMaxDurationInCycles, size_t instrIdx, size_t startCycle, int slot, const std::string &instrumentName)
{
	comment(QL_SS2S(
		"  # slot=" << slot
		<< ", instrument='" << instrumentName << "'"
		<< ": lastEndCycle=" << lastEndCycle[instrIdx]
		<< ", startCycle=" << startCycle
		<< ", instrMaxDurationInCycles=" << instrMaxDurationInCycles
	));

	padToCycle(instrIdx, startCycle, slot, instrumentName);

	// emit code for slot output
	if(condGateMap.empty()) {	// all groups unconditional
		emit(
			slot,
			"seq_out",
			QL_SS2S("0x" << std::hex << std::setfill('0') << std::setw(8) << digOut << std::dec << "," << instrMaxDurationInCycles),
			QL_SS2S("# cycle " << startCycle << "-" << startCycle + instrMaxDurationInCycles << ": code word/mask on '" << instrumentName + "'")
		);
	} else {	// at least one group conditional
		// configure datapath PL
		int smAddr = 0;		// FIXME:
		int pl = dp.getOrAssignPl(instrIdx, condGateMap);

		dp.emitPl(pl, smAddr, condGateMap, instrIdx, slot);

		// emit code for conditional gate
		emit(
			slot,
			"seq_out_sm",
			QL_SS2S("S" << smAddr << "," << pl << "," << instrMaxDurationInCycles),
			QL_SS2S("# cycle " << startCycle << "-" << startCycle + instrMaxDurationInCycles << ": conditional code word/mask on '" << instrumentName << "'")
		);
	}

	// update lastEndCycle
	lastEndCycle[instrIdx] = startCycle + instrMaxDurationInCycles;
}

void codegen_cc::emitPragma(const Json *pragma, int pragmaSmBit, size_t instrIdx, size_t startCycle, int slot, const std::string &instrumentName)
{
	if(startCycle > lastEndCycle[instrIdx]) {	// i.e. if(!instrHasOutput)
		padToCycle(instrIdx, startCycle, slot, instrumentName);
	}

	// FIXME: the only pragma possible is "break" for now
	int pragmaBreakVal = json_get<int>(*pragma, "break", "pragma of unknown instruction");		// FIXME we don't know which instruction we're dealing with, so better move
	int smAddr = pragmaSmBit/32;	// 'seq_cl_sm' is addressable in 32 bit words
	unsigned int mask = 1 << (pragmaSmBit%32);
	std::string label = pragmaForLabel+"_end";		// FIXME: must match label set in forEnd(), assumes we are actually inside a for loop

	// emit code for pragma "break". NB: code is identical for all instruments
/*
	seq_cl_sm   S<address>          ; pass 32 bit SM-data to Q1 (address depends on mapping of variable c) ...
	move_sm     R0                  ; ... and move to register
	and         R0,<mask>,R1        ; mask also depends on mapping of c
	nop								; register dependency R1
	jlt         R1,1,@loop
*/
	emit(slot, "seq_cl_sm", QL_SS2S("S" << smAddr), QL_SS2S("# 'break if " << pragmaBreakVal << "' on '" << instrumentName << "'"));
	emit(slot, "move_sm", "R0", "");
	emit(slot, "and", QL_SS2S("R0," << mask << "," << "R1"), "");	// results in '0' for 'bit==0' and 'mask' for 'bit==1'
	emit(slot, "nop", "", "");
	if(pragmaBreakVal==0) {
		emit(slot, "jlt", QL_SS2S("R1,1,@" << label), "");
	} else {
		emit(slot, "jgt", QL_SS2S("R1,0,@" << label), "");
	}
}


void codegen_cc::padToCycle(size_t instrIdx, size_t startCycle, int slot, const std::string &instrumentName)
{
    // compute prePadding: time to bridge to align timing
    int prePadding = startCycle - lastEndCycle[instrIdx];
    if(prePadding < 0) {
        QL_EOUT("Inconsistency detected in bundle contents: printing code generated so far");
        showCodeSoFar();
        QL_FATAL(
        	"Inconsistency detected in bundle contents: time travel not yet possible in this version: prePadding=" << prePadding
        	<< ", startCycle=" << startCycle
        	<< ", lastEndCycle=" << lastEndCycle[instrIdx]
        	<< ", instrumentName='" << instrumentName << "'"
        	<< ", instrIdx=" << instrIdx
		);
    }

    if(prePadding > 0) {     // we need to align
        emit(
        	slot,
            "seq_wait",
            QL_SS2S(prePadding),
            QL_SS2S("# cycle " << lastEndCycle[instrIdx] << "-" << startCycle << ": padding on '" << instrumentName+"'")
		);
    }

    // update lastEndCycle
    lastEndCycle[instrIdx] = startCycle;
}


// compute signalValueString, and some meta information, for sd[s] (i.e. one of the signals in the JSON definition of an instruction)
codegen_cc::tCalcSignalValue codegen_cc::calcSignalValue(const settings_cc::tSignalDef &sd, size_t s, const Vec<UInt> &operands, const std::string &iname)
{   tCalcSignalValue ret;
    std::string signalSPath = QL_SS2S(sd.path<<"["<<s<<"]");                   // for JSON error reporting

	/************************************************************************\
	| get signal properties, mapping operand index to qubit
	\************************************************************************/

    // get the operand index & qubit to work on
    ret.operandIdx = json_get<unsigned int>(sd.signal[s], "operand_idx", signalSPath);
    if(ret.operandIdx >= operands.size()) {
        QL_JSON_FATAL(
        	"instruction '" << iname
        	<< "': JSON file defines operand_idx " << ret.operandIdx
        	<< ", but only " << operands.size()
        	<< " operands were provided (correct JSON, or provide enough operands)"
		); // FIXME: add offending statement
    }
    UInt qubit = operands[ret.operandIdx];

    // get signal value
    const Json instructionSignalValue = json_get<const Json>(sd.signal[s], "value", signalSPath);   // NB: json_get<const Json&> unavailable
    std::string sv = QL_SS2S(instructionSignalValue);   // serialize/stream instructionSignalValue into std::string

	// get instruction signal type (e.g. "mw", "flux", etc)
    // NB: instructionSignalType is different from "instruction/type" provided by find_instruction_type, although some identical strings are used). NB: that key is no longer used by the 'core' of OpenQL
    std::string instructionSignalType = json_get<std::string>(sd.signal[s], "type", signalSPath);

	/************************************************************************\
	| map signal type for qubit to instrument & group
	\************************************************************************/

    // find signalInfo, i.e. perform the mapping
    ret.si = settings.findSignalInfoForQubit(instructionSignalType, qubit);

	if(instructionSignalValue.size() == 0) {	// allow empty signal
		ret.signalValueString = "";
	} else {
		// verify signal dimensions
		size_t channelsPergroup = ret.si.ic.controlModeGroupSize;
		if(instructionSignalValue.size() != channelsPergroup) {
			QL_JSON_FATAL(
				"signal dimension mismatch on instruction '" << iname
				<< "' : control mode '" << ret.si.ic.refControlMode
				<< "' requires " <<  channelsPergroup
				<< " signals, but signal '" << signalSPath+"/value"
				<< "' provides " << instructionSignalValue.size()
			);
		}

		// expand macros
		sv = replace_all(sv, "\"", "");   // get rid of quotes
		sv = replace_all(sv, "{gateName}", iname);
		sv = replace_all(sv, "{instrumentName}", ret.si.ic.ii.instrumentName);
		sv = replace_all(sv, "{instrumentGroup}", std::to_string(ret.si.group));
		// FIXME: allow using all qubits involved (in same signalType?, or refer to signal: qubitOfSignal[n]), e.g. qubit[0], qubit[1], qubit[2]
		sv = replace_all(sv, "{qubit}", std::to_string(qubit));
		ret.signalValueString = sv;

		// FIXME: note that the actual contents of the signalValue only become important when we'll do automatic codeword assignment and provide codewordTable to downstream software to assign waveforms to the codewords
	}

    comment(QL_SS2S(
    	"  # slot=" << ret.si.ic.ii.slot
		<< ", instrument='" << ret.si.ic.ii.instrumentName << "'"
		<< ", group=" << ret.si.group
		<< "': signalValue='" << ret.signalValueString << "'"
	));

    return ret;
}


#if !OPT_SUPPORT_STATIC_CODEWORDS
uint32_t codegen_cc::assignCodeword(const std::string &instrumentName, int instrIdx, int group)
{
    uint32_t codeword;
    std::string signalValue = bi->signalValue;

    if(QL_JSON_EXISTS(codewordTable, instrumentName) &&                    	// instrument exists
                    codewordTable[instrumentName].size() > group) {     	// group exists
        bool cwFound = false;
        // try to find signalValue
        Json &myCodewordArray = codewordTable[instrumentName][group];
        for(codeword=0; codeword<myCodewordArray.size() && !cwFound; codeword++) {   // NB: JSON find() doesn't work for arrays
            if(myCodewordArray[codeword] == signalValue) {
                QL_DOUT("signal value found at cw=" << codeword);
                cwFound = true;
            }
        }
        if(!cwFound) {
            std::string msg = QL_SS2S("signal value '" << signalValue
                    << "' not found in group " << group
                    << ", which contains " << myCodewordArray);
            if(mapPreloaded) {
                QL_FATAL("mismatch between preloaded 'backend_cc_map_input_file' and program requirements:" << msg)
            } else {
                QL_DOUT(msg);
                // NB: codeword already contains last used value + 1
                // FIXME: check that number is available
                myCodewordArray[codeword] = signalValue;                    // NB: structure created on demand
            }
        }
    } else {    // new instrument or group
        if(mapPreloaded) {
            QL_FATAL("mismatch between preloaded 'backend_cc_map_input_file' and program requirements: instrument '"
                  << instrumentName << "', group "
                  << group
                  << " not present in file");
        } else {
            codeword = 1;
            codewordTable[instrumentName][group][0] = "";                   // code word 0 is empty
            codewordTable[instrumentName][group][codeword] = signalValue;   // NB: structure created on demand
        }
    }
    return codeword;
}
#endif

} // namespace ql
