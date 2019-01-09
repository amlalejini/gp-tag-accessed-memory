#ifndef PROGRAMMING_SYNTHESIS_EXP_H
#define PROGRAMMING_SYNTHESIS_EXP_H

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <unordered_set>
#include <tuple>

#include "base/Ptr.h"
#include "base/vector.h"
#include "control/Signal.h"
#include "Evolve/World.h"
#include "Evolve/World_select.h"
#include "tools/Random.h"
#include "tools/random_utils.h"
#include "tools/math.h"
#include "tools/string_utils.h"
#include "tools/stats.h"
#include "tools/tuple_utils.h"

#include "TagLinearGP.h"
#include "TagLinearGP_InstLib.h"
#include "TagLinearGP_Utilities.h"

#include "TestCaseSet.h"
#include "Selection.h"
#include "Mutators.h"
#include "ProgOrg.h"

#include "ProgSynthConfig.h"
#include "ProgSynthBenchmarks_InputReps.h"

//////////////////////////////////////////
// --- Todo ---
// - [ ] Implement capacity for 'no search' tag argument treatment
//////////////////////////////////////////

//////////////////////////////////////////
// --- Notes ---
// - This experiment setup is largely modeled/copied from the programming synthesis
//   experiment code in my antagonistic-lexicase github repository.
//////////////////////////////////////////

//////////////////////////////////////////
// - Numeric-related instructions -
//   Add
//   Sub
//   Mult
//   Div
//   Mod
//   TestNumEqu
//   TestNumNEqu
//   TestNumLess
//   TestNumLessTEqu
//   TestNumGreater
//   TestNumGreaterTEqu
//   Floor
//   Not
//   Inc
//   Dec
// - Memory-related instructions -
//   CopyMem
//   SwapMem
//   Input
//   Output
//   CommitGlobal
//   PullGlobal
// - Vector-related instructions -
//   MakeVector
//   VecGet
//   VecSet
//   VecLen
//   VecAppend
//   VecPop
//   VecRemove
//   VecReplaceAll
//   VecIndexOf
//   VecOccurrencesOf
//   VecReverse
//   VecSwapIfLess
//   VecGetFront
//   VecGetBack
//   Foreach
// - String-related instructions -
//   StrLength
//   StrConcat
// - Type-related instructions -
//   IsNum
//   IsStr
//   IsVec
//   TestMemEqu
//   TestMemNEqu
// - Non-module flow control instructions -
//   If
//   IfNot
//   While
//   Countdown
//   Close
//   Break
// - Module-related flow control instructions -
//   Call
//   Routine
//   Return
//   ModuleDef
//////////////////////////////////////////


constexpr size_t TAG_WIDTH = 16;
constexpr size_t MEM_SIZE = TAG_WIDTH;

enum PROGRAM_ARGUMENT_MODE_TYPE { TAG_ONLY=0, NUMERIC_ONLY=1, BOTH=2 };

enum PROBLEM_ID { NumberIO=0,
                  SmallOrLarge,
                  ForLoopIndex,
                  CompareStringLengths,
                  DoubleLetters,
                  CollatzNumbers,
                  ReplaceSpaceWithNewline,
                  StringDifferences,
                  EvenSquares,
                  WallisPi,
                  StringLengthsBackwards,
                  LastIndexOfZero,
                  VectorAverage,
                  CountOdds,
                  MirrorImage,
                  SuperAnagrams,
                  SumOfSquares,
                  VectorsSummed,
                  XWordLines,
                  PigLatin,
                  NegativeToZero,
                  ScrabbleScore,
                  Checksum,
                  Digits,
                  Grade,
                  Median,
                  Smallest,
                  Syllables
};

// Structure to track problem logistics.
struct ProblemInfo {
  PROBLEM_ID id; 
  std::string training_fname;
  std::string testing_fname;
  
  ProblemInfo(PROBLEM_ID _id, const std::string & _training_fname, const std::string & _testing_fname) 
    : id(_id), training_fname(_training_fname), testing_fname(_testing_fname)
  { ; }
  
  ProblemInfo(const ProblemInfo &) = default;
  ProblemInfo(ProblemInfo &&) = default;

  ProblemInfo & operator=(const ProblemInfo &) = default;
  ProblemInfo & operator=(ProblemInfo &&) = default;

  const std::string & GetTestingSetFilename() const { return testing_fname; }
  const std::string & GetTrainingSetFilename() const { return training_fname; }

};

// Map of problem name to problem information.
std::unordered_map<std::string, ProblemInfo> problems = {
  {"number-io", {PROBLEM_ID::NumberIO, "training-examples-number-io.csv", "testing-examples-number-io.csv"}},
  {"small-or-large", {PROBLEM_ID::SmallOrLarge, "training-examples-small-or-large.csv", "testing-examples-small-or-large.csv"}},
  {"for-loop-index", {PROBLEM_ID::ForLoopIndex, "training-examples-for-loop-index.csv", "testing-examples-for-loop-index.csv"}},
  {"compare-string-lengths", {PROBLEM_ID::CompareStringLengths, "training-examples-compare-string-lengths.csv", "testing-examples-compare-string-lengths.csv"}},
  {"collatz-numbers", {PROBLEM_ID::CollatzNumbers, "training-examples-collatz-numbers.csv", "testing-examples-collatz-numbers.csv"}},
  {"string-lengths-backwards", {PROBLEM_ID::StringLengthsBackwards, "training-examples-string-lengths-backwards.csv", "testing-examples-string-lengths-backwards.csv"}},
  {"last-index-of-zero", {PROBLEM_ID::LastIndexOfZero, "training-examples-last-index-of-zero.csv", "testing-examples-last-index-of-zero.csv"}},
  {"count-odds", {PROBLEM_ID::CountOdds, "training-examples-count-odds.csv", "testing-examples-count-odds.csv"}},
  {"mirror-image", {PROBLEM_ID::MirrorImage, "training-examples-mirror-image.csv", "testing-examples-mirror-image.csv"}},
  {"vectors-summed", {PROBLEM_ID::VectorsSummed, "training-examples-vectors-summed.csv", "testing-examples-vectors-summed.csv"}},
  {"sum-of-squares", {PROBLEM_ID::SumOfSquares, "training-examples-sum-of-squares.csv", "testing-examples-sum-of-squares.csv"}},
  {"vector-average", {PROBLEM_ID::VectorAverage, "training-examples-vector-average.csv", "testing-examples-vector-average.csv"}},
  {"median", {PROBLEM_ID::Median, "training-examples-median.csv", "testing-examples-median.csv"}},
  {"smallest", {PROBLEM_ID::Smallest, "training-examples-smallest.csv", "testing-examples-smallest.csv"}},
  {"grade", {PROBLEM_ID::Grade, "training-examples-grade.csv", "testing-examples-grade.csv"}}
};

class ProgramSynthesisExperiment {
public:

  using hardware_t = typename TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using inst_lib_t = typename TagLGP::InstLib<hardware_t>;
  using inst_t = typename hardware_t::inst_t;

  using prog_org_t = ProgOrg<TAG_WIDTH>;
  using prog_org_phen_t = typename prog_org_t::Phenotype;
  using prog_org_gen_t = typename prog_org_t::genome_t;
  using prog_taxon_t = typename emp::Systematics<prog_org_t, prog_org_gen_t>::taxon_t;

  using prog_world_t = emp::World<prog_org_t>;

private:

  // Configuration variables
  int SEED;
  size_t GENERATIONS;
  size_t PROGRAM_ARGUMENT_MODE;
  size_t PROG_POP_SIZE;
  std::string PROBLEM;
  std::string BENCHMARK_DATA_DIR;

  size_t MIN_PROG_SIZE;
  size_t MAX_PROG_SIZE;
  size_t PROG_EVAL_TIME;
  double PROG_MUT__PER_BIT_FLIP;
  double PROG_MUT__PER_NUMERIC_ARG_SUB;
  double PROG_MUT__PER_INST_SUB;
  double PROG_MUT__PER_INST_INS;
  double PROG_MUT__PER_INST_DEL;
  double PROG_MUT__PER_PROG_SLIP;
  double PROG_MUT__PER_MOD_DUP;
  double PROG_MUT__PER_MOD_DEL;

  double MIN_TAG_SPECIFICITY;
  size_t MAX_CALL_DEPTH;

  // Experiment variables
  bool setup;
  size_t update;

  size_t dominant_prog_id;

  emp::Ptr<emp::Random> random;

  emp::Ptr<inst_lib_t> inst_lib;
  emp::Ptr<hardware_t> eval_hardware;

  size_t eval_time;

  size_t smallest_prog_solution_size;
  bool solution_found;
  size_t update_first_solution_found;

  emp::BitSet<TAG_WIDTH> call_tag;

  emp::Ptr<prog_world_t> prog_world;
  // emp::Ptr<emp::Systematics<prog_org_t, prog_org_gen_t>> prog_genotypic_systematics;

  emp::vector<std::function<double(prog_org_t &)>> lexicase_prog_fit_set; /// Fitness function set for performing lexicase selection on programs.

  // Problem utilities
  ProblemUtilities_NumberIO prob_utils_NumberIO;
  ProblemUtilities_SmallOrLarge prob_utils_SmallOrLarge;
  ProblemUtilities_ForLoopIndex prob_utils_ForLoopIndex;
  ProblemUtilities_CompareStringLengths prob_utils_CompareStringLengths;
  ProblemUtilities_DoubleLetters prob_utils_DoubleLetters;
  ProblemUtilities_CollatzNumbers prob_utils_CollatzNumbers;
  ProblemUtilities_ReplaceSpaceWithNewline prob_utils_ReplaceSpaceWithNewline;
  ProblemUtilities_StringDifferences prob_utils_StringDifferences;
  ProblemUtilities_EvenSquares prob_utils_EvenSquares;
  ProblemUtilities_WallisPi prob_utils_WallisPi;
  ProblemUtilities_StringLengthsBackwards prob_utils_StringLengthsBackwards;
  ProblemUtilities_LastIndexOfZero prob_utils_LastIndexOfZero;
  ProblemUtilities_VectorAverage prob_utils_VectorAverage;
  ProblemUtilities_CountOdds prob_utils_CountOdds;
  ProblemUtilities_MirrorImage prob_utils_MirrorImage;
  ProblemUtilities_SuperAnagrams prob_utils_SuperAnagrams;
  ProblemUtilities_SumOfSquares prob_utils_SumOfSquares;
  ProblemUtilities_VectorsSummed prob_utils_VectorsSummed;
  ProblemUtilities_XWordLines prob_utils_XWordLines;
  ProblemUtilities_PigLatin prob_utils_PigLatin;
  ProblemUtilities_NegativeToZero prob_utils_NegativeToZero;
  ProblemUtilities_ScrabbleScore prob_utils_ScrabbleScore;
  ProblemUtilities_Checksum prob_utils_Checksum;
  ProblemUtilities_Digits prob_utils_Digits;
  ProblemUtilities_Grade prob_utils_Grade;
  ProblemUtilities_Median prob_utils_Median;
  ProblemUtilities_Smallest prob_utils_Smallest;
  ProblemUtilities_Syllables prob_utils_Syllables;

  emp::Ptr<emp::DataFile> solution_file;
  
  /// Struct to track evaluation test result information.
  struct TestResult {
    double score;   ///< Program score on test.
    bool pass;      ///< Did program pass the test?
    bool sub;       ///< Did program submit output for the test?
    TestResult(double sc=0, bool p=false, bool sb=false) : score(sc), pass(p), sub(sb) { ; }
  };

  /// Utility for managing program and test information during evaluations.
  /// - Especially useful during population snapshots, solution validations, etc.
  struct EvalUtil {
    size_t current_programID;
    size_t current_testID;

    EvalUtil(size_t pID=0, size_t tID=0) : current_programID(pID), current_testID(tID) { ; }
  } eval_util;

  // Experiment signals
  emp::Signal<void(void)> do_evaluation_sig;
  emp::Signal<void(void)> do_selection_sig;
  emp::Signal<void(void)> do_update_sig;  

  emp::Signal<void(void)> end_setup_sig;

  // Program evaluation signals
  emp::Signal<void(prog_org_t &)> begin_program_eval;   /// Begin evaluating a program on a test case set.
  emp::Signal<void(prog_org_t &)> end_program_eval;     /// Finish evaluating a program on a test case set.

  emp::Signal<void(prog_org_t &)> begin_program_test; /// Begin evaluating a program on an individual test.
  emp::Signal<void(prog_org_t &)> do_program_test;    /// Do single-test evaluation with given program on given test.
  emp::Signal<void(prog_org_t &)> end_program_test;   /// Finish evaluating a program on an individual test.

  emp::Signal<void(prog_org_t &)> do_program_advance; /// Advance virtual hardware by one time step.

  // std::function<void(prog_org_t &)> ValidateProgram

  // Internal functions
  void InitConfigs(const ProgramSynthesisConfig & config);
  
  void InitProgPop_Random();

  void AddDefaultInstructions_TagArgs();
  void AddDefaultInstructions_TagArgs_NoTypeSearch();
  void AddDefaultInstructions_NumArgs();

  void AddVectorInstructions_TagArgs();
  void AddVectorInstructions_TagArgs_NoTypeSearch();
  void AddVectorInstructions_NumArgs();

  void AddStringInstructions_TagArgs();
  void AddStringInstructions_TagArgs_NoTypeSearch();
  void AddStringInstructions_NumArgs();
  
  void SetupHardware();
  void SetupEvaluation();
  void SetupSelection();
  void SetupMutation();
  void SetupFitFuns();

  void SetupProblem();
  void SetupProblem_NumberIO();
  void SetupProblem_SmallOrLarge();
  void SetupProblem_ForLoopIndex();
  void SetupProblem_CompareStringLengths();
  void SetupProblem_DoubleLetters();
  void SetupProblem_CollatzNumbers();
  void SetupProblem_ReplaceSpaceWithNewline();
  void SetupProblem_StringDifferences();
  void SetupProblem_EvenSquares();
  void SetupProblem_WallisPi();
  void SetupProblem_StringLengthsBackwards();
  void SetupProblem_LastIndexOfZero();
  void SetupProblem_VectorAverage();
  void SetupProblem_CountOdds();
  void SetupProblem_MirrorImage();
  void SetupProblem_SuperAnagrams();
  void SetupProblem_SumOfSquares();
  void SetupProblem_VectorsSummed();
  void SetupProblem_XWordLines();
  void SetupProblem_PigLatin();
  void SetupProblem_NegativeToZero();
  void SetupProblem_ScrabbleScore();
  void SetupProblem_Checksum();
  void SetupProblem_Digits();
  void SetupProblem_Grade();
  void SetupProblem_Median();
  void SetupProblem_Smallest();
  void SetupProblem_Syllables();

public:

  ProgramSynthesisExperiment()
    : setup(false), update(0), solution_found(false) 
  {
    std::cout << "Problem information:" << std::endl;
    for (const auto & info : problems) {
      std::cout << "  - Problem name: " << info.first << std::endl;
      std::cout << "    - Training examples file: " << info.second.GetTrainingSetFilename() << std::endl;
      std::cout << "    - Testing examples file: " << info.second.GetTestingSetFilename() << std::endl;
    }
  }

  ~ProgramSynthesisExperiment() {
    if (setup) {
      // todo
    }
  }

  /// Configure the experiment.
  void Setup(const ProgramSynthesisConfig & config);

  /// Run the experiment start->finish.
  void Run();

  /// Progress the experiment by a single time step (generation).
  /// (1) evaluate the population
  /// (2) select individuals for reproduction
  /// (3) update world(s)
  void RunStep();

};

/// Configure the experiment.
void ProgramSynthesisExperiment::Setup(const ProgramSynthesisConfig & config) {
  std::cout << "Running Program synthesis experiment setup." << std::endl;
  emp_assert(setup == false, "Can only run setup once because, lazy.");
  // Localize experiment configuration.
  InitConfigs(config);

  if (setup) {
    // Fail!
    std::cout << "Setup is only allowed once per experiment! Exiting." << std::endl;
    exit(-1); 
  }

  // Allocate memory for random number generator and program world.
  random = emp::NewPtr<emp::Random>(SEED);
  prog_world = emp::NewPtr<prog_world_t>(*random, "Program World");

  // Initialize solution tracking variables.
  smallest_prog_solution_size = MAX_PROG_SIZE + 1;
  solution_found = false;
  update_first_solution_found = GENERATIONS + 1;
  dominant_prog_id = 0;

  // Configure the program world.
  prog_world->SetPopStruct_Mixed(true);

  // Configure how program population should be initialized.
  end_setup_sig.AddAction([this]() {
    // Initialize the program population.
    InitProgPop_Random();
    std::cout << "Done initializing program population. Program population size = " << prog_world->GetSize() << std::endl;
  });

  // Configure the on update signal.
  do_update_sig.AddAction([this]() {
    std::cout << "Update: " << update << "; ";
    std::cout << "best program score: " << prog_world->CalcFitnessID(dominant_prog_id) << "; ";
    std::cout << "solution found? " << solution_found << "; ";
    std::cout << "smallest solution? " << smallest_prog_solution_size << std::endl;

    // if (update % SNAPSHOT_INTERVAL == 0 || update == GENERATIONS) do_pop_snapshot_sig.Trigger(); 

    prog_world->Update();
    prog_world->ClearCache();
  });

  // Setup the virtual hardware
  std::cout << "==== EXPERIMENT SETUP => evaluation hardware ====" << std::endl;
  SetupHardware(); 

  // Setup problem that we're evolving programs to solve. The particular problem
  // we setup depends on experiment configuration.
  std::cout << "==== EXPERIMENT SETUP => problem ====" << std::endl;
  // SetupProblem(); -- todo --


  // Trigger end of setup.
  end_setup_sig.Trigger();
}

void ProgramSynthesisExperiment::SetupHardware() {
  // Create a new instruction library.
  inst_lib = emp::NewPtr<inst_lib_t>();
  // Create evaluation hardware.
  eval_hardware = emp::NewPtr<hardware_t>(inst_lib, random);
  // Configure the CPU.
  eval_hardware->SetMemSize(MEM_SIZE);
  eval_hardware->SetMinTagSpecificity(MIN_TAG_SPECIFICITY);  // Configure minimum tag specificity required for tag-based referencing.
  eval_hardware->SetMaxCallDepth(MAX_CALL_DEPTH);            // Configure maximum depth of call stack (recursion limit).
  eval_hardware->SetMemTags(GenHadamardMatrix<TAG_WIDTH>()); // Configure memory location tags. Use Hadamard matrix for given TAG_WIDTH.

  // Configure call tag (tag used to call initial module during test evaluation).
  call_tag.Clear(); // Set initial call tag to all 0s.

  // What do we do at the beginning of program evaluation?
  begin_program_eval.AddAction([this](prog_org_t & prog_org) {
    eval_hardware->Reset();
    eval_hardware->SetProgram(prog_org.GetGenome());
  });

  // What should we do to the hardware after program evaluation?
  // - Currently, nothing.

  // What do we do before running a program on a single test?
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    eval_hardware->ResetHardware();
    eval_hardware->CallModule(call_tag, MIN_TAG_SPECIFICITY, true, false); /// If we're not using any modules (no module instructions), this will always call the default module!
  });

  // How do we 'do' a program test?
  // - For a specified evaluation time, advance the evaluation hardware.
  do_program_test.AddAction([this](prog_org_t & prog_org) {
    // std::cout << "--- DO PROGRAM TEST ---" << std::endl;
    // std::cout << "==== Initial hardware state ====" << std::endl;
    // eval_hardware->PrintHardwareState();
    for (eval_time = 0; eval_time < PROG_EVAL_TIME; ++eval_time) {
      // std::cout << "==== Time = " << eval_time << "==== " << std::endl;
      do_program_advance.Trigger(prog_org);
      // eval_hardware->PrintHardwareState();
      if (eval_hardware->GetCallStackSize() == 0) break; // If call stack is ever completely empty, program is done early.
    }
    // exit(-1);
  });

  // How do we advance the evaluation hardware?
  do_program_advance.AddAction([this](prog_org_t &) {
    eval_hardware->SingleProcess();
  });

  // Setup default instruction set.
  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      AddDefaultInstructions_TagArgs();
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      AddDefaultInstructions_NumArgs();
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      AddDefaultInstructions_NumArgs();
      AddDefaultInstructions_TagArgs();
      break;
    }
    default: {
      std::cout << "Unrecognized PROGRAM_ARGUMENT_MODE (" << PROGRAM_ARGUMENT_MODE << "). Exiting." << std::endl;
      exit(-1); 
    }
  }

}

void ProgramSynthesisExperiment::SetupProblem() {

}

void ProgramSynthesisExperiment::SetupEvaluation() {

}

void ProgramSynthesisExperiment::SetupSelection() {

}

void ProgramSynthesisExperiment::SetupMutation() {

}

void ProgramSynthesisExperiment::SetupFitFuns() {

}


void ProgramSynthesisExperiment::InitConfigs(const ProgramSynthesisConfig & config) {

  SEED = config.SEED();
  GENERATIONS = config.GENERATIONS();
  PROGRAM_ARGUMENT_MODE = config.PROGRAM_ARGUMENT_MODE();
  PROG_POP_SIZE = config.PROG_POP_SIZE();
  PROBLEM = config.PROBLEM();
  BENCHMARK_DATA_DIR = config.BENCHMARK_DATA_DIR();

  MIN_PROG_SIZE = config.MIN_PROG_SIZE();
  MAX_PROG_SIZE = config.MAX_PROG_SIZE();
  PROG_EVAL_TIME = config.PROG_EVAL_TIME();
  PROG_MUT__PER_BIT_FLIP = config.PROG_MUT__PER_BIT_FLIP();
  PROG_MUT__PER_NUMERIC_ARG_SUB = config.PROG_MUT__PER_NUMERIC_ARG_SUB();
  PROG_MUT__PER_INST_SUB = config.PROG_MUT__PER_INST_SUB();
  PROG_MUT__PER_INST_INS = config.PROG_MUT__PER_INST_INS();
  PROG_MUT__PER_INST_DEL = config.PROG_MUT__PER_INST_DEL();
  PROG_MUT__PER_PROG_SLIP = config.PROG_MUT__PER_PROG_SLIP();
  PROG_MUT__PER_MOD_DUP = config.PROG_MUT__PER_MOD_DUP();
  PROG_MUT__PER_MOD_DEL = config.PROG_MUT__PER_MOD_DEL();
  
  MIN_TAG_SPECIFICITY = config.MIN_TAG_SPECIFICITY();
  MAX_CALL_DEPTH = config.MAX_CALL_DEPTH();

}

void ProgramSynthesisExperiment::InitProgPop_Random() {
  std::cout << "Randomly initializing program population." << std::endl;
  for (size_t i = 0; i < PROG_POP_SIZE; ++i) {
    switch (PROGRAM_ARGUMENT_MODE) {
      case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
        prog_world->Inject(TagLGP::GenRandTagGPProgram(*random, inst_lib, MIN_PROG_SIZE, MAX_PROG_SIZE), 1);
        break;
      }
      case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
        prog_world->Inject(TagLGP::GenRandTagGPProgram_NumArgs(*random, inst_lib, MEM_SIZE-1, MIN_PROG_SIZE, MAX_PROG_SIZE), 1);
        break;
      }
      case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
        prog_world->Inject(TagLGP::GenRandTagGPProgram_TagAndNumArgs(*random, inst_lib, MEM_SIZE-1, MIN_PROG_SIZE, MAX_PROG_SIZE), 1);
        break;
      }
      default: {
        std::cout << "Unrecognized PROGRAM_ARGUMENT_MODE (" << PROGRAM_ARGUMENT_MODE << "). Exiting." << std::endl;
        exit(-1); 
      }
    }
  }
}

/// Add instructions (expecting tag arguments) that are shared across all problems.
void ProgramSynthesisExperiment::AddDefaultInstructions_TagArgs() {
  // - Numeric-related instructions -
  //   Add
  //   Sub
  //   Mult
  //   Div
  //   Mod
  //   TestNumEqu
  //   TestNumNEqu
  //   TestNumLess
  //   TestNumLessTEqu
  //   TestNumGreater
  //   TestNumGreaterTEqu
  //   Floor
  //   Not
  //   Inc
  //   Dec
  inst_lib->AddInst("Add-Tag", hardware_t::Inst_Add, 3, "wmemANY[C] = wmemNUM[A] + wmemNUM[B]");
  inst_lib->AddInst("Sub-Tag", hardware_t::Inst_Sub, 3, "wmemANY[C] = wmemNUM[A] - wmemNUM[B]");
  inst_lib->AddInst("Mult-Tag", hardware_t::Inst_Mult, 3, "wmemANY[C] = wmemNUM[A] * wmemNUM[B]");
  inst_lib->AddInst("Div-Tag", hardware_t::Inst_Div, 3, "if (wmemNUM[B] != 0) wmemANY[C] = wmemNUM[A] / wmemNUM[B]; else NOP");
  inst_lib->AddInst("Mod-Tag", hardware_t::Inst_Mod, 3, "if (wmemNUM[B] != 0) wmemANY[C] = int(wmemNUM[A]) % int(wmemNUM[B]); else NOP");
  inst_lib->AddInst("TestNumEqu-Tag", hardware_t::Inst_TestNumEqu, 3, "wmemANY[C] = wmemNUM[A] == wmemNUM[B]");
  inst_lib->AddInst("TestNumNEqu-Tag", hardware_t::Inst_TestNumNEqu, 3, "wmemANY[C] = wmemNUM[A] != wmemNUM[B]");
  inst_lib->AddInst("TestNumLess-Tag", hardware_t::Inst_TestNumLess, 3, "wmemANY[C] = wmemNUM[A] < wmemNUM[B]");
  inst_lib->AddInst("TestNumLessTEqu-Tag", hardware_t::Inst_TestNumLessTEqu, 3, "wmemANY[C] = wmemNUM[A] <= wmemNUM[B]");
  inst_lib->AddInst("TestNumGreater-Tag", hardware_t::Inst_TestNumGreater, 3, "wmemANY[C] = wmemNUM[A] > wmemNUM[B]");
  inst_lib->AddInst("TestNumGreaterTEqu-Tag", hardware_t::Inst_TestNumGreaterTEqu, 3, "wmemANY[C] = wmemNUM[A] >= wmemNUM[B]");
  inst_lib->AddInst("Floor-Tag", hardware_t::Inst_Floor, 1, "wmemNUM[A] = floor(wmemNUM[A])");
  inst_lib->AddInst("Not-Tag", hardware_t::Inst_Not, 1, "wmemNUM[A] = !wmemNUM[A]"); 
  inst_lib->AddInst("Inc-Tag", hardware_t::Inst_Inc, 1, "wmemNUM[A] = wmemNUM[A] + 1");
  inst_lib->AddInst("Dec-Tag", hardware_t::Inst_Dec, 1, "wmemNUM[A] = wmemNUM[A] - 1");

  // - Memory-related instructions -
  //   CopyMem
  //   SwapMem
  //   Input
  //   Output
  //   CommitGlobal
  //   PullGlobal
  inst_lib->AddInst("CopyMem-Tag", hardware_t::Inst_CopyMem, 2, "wmemANY[B] = wmemANY[A] // Copy mem[A] to mem[B]");
  inst_lib->AddInst("SwapMem-Tag", hardware_t::Inst_SwapMem, 2, "swap(wmemANY[A], wmemANY[B])");
  inst_lib->AddInst("Input-Tag", hardware_t::Inst_Input, 2, "wmemANY[B] = imemANY[A]");
  inst_lib->AddInst("Output-Tag", hardware_t::Inst_Output, 2, "omemANY[B] = wmemANY[A]");
  inst_lib->AddInst("CommitGlobal-Tag", hardware_t::Inst_CommitGlobal, 2, "gmemANY[B] = wmemANY[A]");
  inst_lib->AddInst("PullGlobal-Tag", hardware_t::Inst_PullGlobal, 2, "wmemANY[B] = gmemANY[A]");
  // inst_lib->AddInst("TestMemEqu-Tag", hardware_t::Inst_TestMemEqu, 3, "wmemANY[C] = wmemANY[A] == wmemANY[B]");
  // inst_lib->AddInst("TestMemNEqu-Tag", hardware_t::Inst_TestMemNEqu, 3, "wmemANY[C] = wmemANY[A] != wmemANY[B]");

  // - Non-module flow control instructions -
  //   If
  //   IfNot
  //   While
  //   Countdown
  //   Close
  //   Break
  inst_lib->AddInst("If-Tag", hardware_t::Inst_If, 1, "Execute next flow if(wmemANY[A]) // To be true, mem loc must be non-zero number", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("IfNot-Tag", hardware_t::Inst_IfNot, 1, "Execute next flow if(!wmemANY[A])", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("While-Tag", hardware_t::Inst_While, 1, "While loop over wmemANY[A]", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Countdown-Tag", hardware_t::Inst_Countdown, 1, "Countdown loop with wmemANY as index.", {inst_lib_t::InstProperty::BEGIN_FLOW});
  // inst_lib->AddInst("Foreach-Tag", hardware_t::Inst_Foreach, 2, "For each thing in wmemVEC[B]", {inst_lib_t::InstProperty::BEGIN_FLOW});
  
  // The below instructions take no arguments, check if instruction in library first.
  if (!inst_lib->IsInst("Close")) {
    inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "Close flow", {inst_lib_t::InstProperty::END_FLOW});
  }
  if (!inst_lib->IsInst("Break")) {
    inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "Break current flow");
  }
}

/// Add instructions (expecting tag arguments) that are shared across all problems.
void ProgramSynthesisExperiment::AddDefaultInstructions_NumArgs() {
  // - Numeric-related instructions -
  //   Add
  //   Sub
  //   Mult
  //   Div
  //   Mod
  //   TestNumEqu
  //   TestNumNEqu
  //   TestNumLess
  //   TestNumLessTEqu
  //   TestNumGreater
  //   TestNumGreaterTEqu
  //   Floor
  //   Not
  //   Inc
  //   Dec
  inst_lib->AddInst("Add-Num", hardware_t::Inst_Add__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] + wmemNUM[B]");
  inst_lib->AddInst("Sub-Num", hardware_t::Inst_Sub__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] - wmemNUM[B]");
  inst_lib->AddInst("Mult-Num", hardware_t::Inst_Mult__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] * wmemNUM[B]");
  inst_lib->AddInst("Div-Num", hardware_t::Inst_Div__NUM_ARGS, 3, "if (wmemNUM[B] != 0) wmemANY[C] = wmemNUM[A] / wmemNUM[B]; else NOP");
  inst_lib->AddInst("Mod-Num", hardware_t::Inst_Mod__NUM_ARGS, 3, "if (wmemNUM[B] != 0) wmemANY[C] = int(wmemNUM[A]) % int(wmemNUM[B]); else NOP");
  inst_lib->AddInst("TestNumEqu-Num", hardware_t::Inst_TestNumEqu__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] == wmemNUM[B]");
  inst_lib->AddInst("TestNumNEqu-Num", hardware_t::Inst_TestNumNEqu__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] != wmemNUM[B]");
  inst_lib->AddInst("TestNumLess-Num", hardware_t::Inst_TestNumLess__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] < wmemNUM[B]");
  inst_lib->AddInst("TestNumLessTEqu-Num", hardware_t::Inst_TestNumLessTEqu__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] <= wmemNUM[B]");
  inst_lib->AddInst("TestNumGreater-Num", hardware_t::Inst_TestNumGreater__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] > wmemNUM[B]");
  inst_lib->AddInst("TestNumGreaterTEqu-Num", hardware_t::Inst_TestNumGreaterTEqu__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] >= wmemNUM[B]");
  inst_lib->AddInst("Floor-Num", hardware_t::Inst_Floor__NUM_ARGS, 1, "wmemNUM[A] = floor(wmemNUM[A])");
  inst_lib->AddInst("Not-Num", hardware_t::Inst_Not__NUM_ARGS, 1, "wmemNUM[A] = !wmemNUM[A]"); 
  inst_lib->AddInst("Inc-Num", hardware_t::Inst_Inc__NUM_ARGS, 1, "wmemNUM[A] = wmemNUM[A] + 1");
  inst_lib->AddInst("Dec-Num", hardware_t::Inst_Dec__NUM_ARGS, 1, "wmemNUM[A] = wmemNUM[A] - 1");

  // - Memory-related instructions -
  //   CopyMem
  //   SwapMem
  //   Input
  //   Output
  //   CommitGlobal
  //   PullGlobal
  inst_lib->AddInst("CopyMem-Num", hardware_t::Inst_CopyMem__NUM_ARGS, 2, "wmemANY[B] = wmemANY[A] // Copy mem[A] to mem[B]");
  inst_lib->AddInst("SwapMem-Num", hardware_t::Inst_SwapMem__NUM_ARGS, 2, "swap(wmemANY[A], wmemANY[B])");
  inst_lib->AddInst("Input-Num", hardware_t::Inst_Input__NUM_ARGS, 2, "wmemANY[B] = imemANY[A]");
  inst_lib->AddInst("Output-Num", hardware_t::Inst_Output__NUM_ARGS, 2, "omemANY[B] = wmemANY[A]");
  inst_lib->AddInst("CommitGlobal-Num", hardware_t::Inst_CommitGlobal__NUM_ARGS, 2, "gmemANY[B] = wmemANY[A]");
  inst_lib->AddInst("PullGlobal-Num", hardware_t::Inst_PullGlobal__NUM_ARGS, 2, "wmemANY[B] = gmemANY[A]");
  // inst_lib->AddInst("TestMemEqu-Num", hardware_t::Inst_TestMemEqu__NUM_ARGS, 3, "wmemANY[C] = wmemANY[A] == wmemANY[B]");
  // inst_lib->AddInst("TestMemNEqu-Num", hardware_t::Inst_TestMemNEqu__NUM_ARGS, 3, "wmemANY[C] = wmemANY[A] != wmemANY[B]");

  // - Non-module flow control instructions -
  //   If
  //   IfNot
  //   While
  //   Countdown
  //   Close
  //   Break
  inst_lib->AddInst("If-Num", hardware_t::Inst_If__NUM_ARGS, 1, "Execute next flow if(wmemANY[A]) // To be true, mem loc must be non-zero number", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("IfNot-Num", hardware_t::Inst_IfNot__NUM_ARGS, 1, "Execute next flow if(!wmemANY[A])", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("While-Num", hardware_t::Inst_While__NUM_ARGS, 1, "While loop over wmemANY[A]", {inst_lib_t::InstProperty::BEGIN_FLOW});
  inst_lib->AddInst("Countdown-Num", hardware_t::Inst_Countdown__NUM_ARGS, 1, "Countdown loop with wmemANY as index.", {inst_lib_t::InstProperty::BEGIN_FLOW});
  // inst_lib->AddInst("Foreach-Num", hardware_t::Inst_Foreach__NUM_ARGS, 2, "For each thing in wmemVEC[B]", {inst_lib_t::InstProperty::BEGIN_FLOW});
  
  // The below instructions take no arguments, check if instruction in library first.
  if (!inst_lib->IsInst("Close")) {
    inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "Close flow", {inst_lib_t::InstProperty::END_FLOW});
  }
  if (!inst_lib->IsInst("Break")) {
    inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "Break current flow");
  }
}

#endif