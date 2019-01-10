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

// Some useful constants for each problem.
constexpr double PROB_NUMBER_IO__DOUBLE_MIN = -100.0;
constexpr double PROB_NUMBER_IO__DOUBLE_MAX =  100.0;
constexpr int PROB_NUMBER_IO__INT_MIN = -100;
constexpr int PROB_NUMBER_IO__INT_MAX =  100;

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
  {"number-io", {PROBLEM_ID::NumberIO, "training-examples-number-io.csv", "testing-examples-number-io.csv"}}
  // {"small-or-large", {PROBLEM_ID::SmallOrLarge, "training-examples-small-or-large.csv", "testing-examples-small-or-large.csv"}},
  // {"for-loop-index", {PROBLEM_ID::ForLoopIndex, "training-examples-for-loop-index.csv", "testing-examples-for-loop-index.csv"}},
  // {"compare-string-lengths", {PROBLEM_ID::CompareStringLengths, "training-examples-compare-string-lengths.csv", "testing-examples-compare-string-lengths.csv"}},
  // {"collatz-numbers", {PROBLEM_ID::CollatzNumbers, "training-examples-collatz-numbers.csv", "testing-examples-collatz-numbers.csv"}},
  // {"string-lengths-backwards", {PROBLEM_ID::StringLengthsBackwards, "training-examples-string-lengths-backwards.csv", "testing-examples-string-lengths-backwards.csv"}},
  // {"last-index-of-zero", {PROBLEM_ID::LastIndexOfZero, "training-examples-last-index-of-zero.csv", "testing-examples-last-index-of-zero.csv"}},
  // {"count-odds", {PROBLEM_ID::CountOdds, "training-examples-count-odds.csv", "testing-examples-count-odds.csv"}},
  // {"mirror-image", {PROBLEM_ID::MirrorImage, "training-examples-mirror-image.csv", "testing-examples-mirror-image.csv"}},
  // {"vectors-summed", {PROBLEM_ID::VectorsSummed, "training-examples-vectors-summed.csv", "testing-examples-vectors-summed.csv"}},
  // {"sum-of-squares", {PROBLEM_ID::SumOfSquares, "training-examples-sum-of-squares.csv", "testing-examples-sum-of-squares.csv"}},
  // {"vector-average", {PROBLEM_ID::VectorAverage, "training-examples-vector-average.csv", "testing-examples-vector-average.csv"}},
  // {"median", {PROBLEM_ID::Median, "training-examples-median.csv", "testing-examples-median.csv"}},
  // {"smallest", {PROBLEM_ID::Smallest, "training-examples-smallest.csv", "testing-examples-smallest.csv"}},
  // {"grade", {PROBLEM_ID::Grade, "training-examples-grade.csv", "testing-examples-grade.csv"}}
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

  bool USE_MODULES;
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

  size_t LEXICASE_MAX_FUNS;

  // Experiment variables
  size_t TRAINING_SET_SIZE;
  size_t TESTING_SET_SIZE;

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
  TagLGPMutator<TAG_WIDTH> prog_mutator;

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
    bool use_training_set;

    EvalUtil(size_t pID=0, size_t tID=0) : current_programID(pID), current_testID(tID), use_training_set(true) { ; }
  } eval_util;

  // Experiment signals
  emp::Signal<void(void)> do_evaluation_sig;
  emp::Signal<void(void)> do_selection_sig;
  emp::Signal<void(void)> do_update_sig;  

  emp::Signal<void(void)> end_setup_sig;

  // Program evaluation signals
  emp::Signal<void(prog_org_t &)> begin_program_eval;   ///< Begin evaluating a program on a test case set.
  emp::Signal<void(prog_org_t &)> end_program_eval;     ///< Finish evaluating a program on a test case set.

  emp::Signal<void(prog_org_t &)> begin_program_test;   ///< Begin evaluating a program on an individual test.
  emp::Signal<void(prog_org_t &)> do_program_test;      ///< Do single-test evaluation with given program on given test.
  emp::Signal<void(prog_org_t &)> end_program_test;     ///< Finish evaluating a program on an individual test.

  emp::Signal<void(prog_org_t &)> do_program_advance;   ///< Advance virtual hardware by one time step.

  // std::function<void(prog_org_t &)> ValidateProgram

  std::function<TestResult(prog_org_t &)> CalcProgramResultOnTest; ///< Evaluate given program on given training case.
  
  // Internal functions
  void InitConfigs(const ProgramSynthesisConfig & config);
  
  void InitProgPop_Random();

  void AddDefaultInstructions();
  void AddDefaultInstructions_TagArgs();
  void AddDefaultInstructions_TagArgs_NoTypeSearch();
  void AddDefaultInstructions_NumArgs();

  void AddVectorInstructions();
  void AddVectorInstructions_TagArgs();
  void AddVectorInstructions_TagArgs_NoTypeSearch();
  void AddVectorInstructions_NumArgs();

  void AddStringInstructions();
  void AddStringInstructions_TagArgs();
  void AddStringInstructions_TagArgs_NoTypeSearch();
  void AddStringInstructions_NumArgs();

  void AddNumericTerminals(size_t min, size_t max);
  void AddNumericTerminals_TagArgs(size_t min, size_t max);
  void AddNumericTerminals_TagArgs_NoTypeSearch(size_t min, size_t max);
  void AddNumericTerminals_NumArgs(size_t min, size_t max);
  
  void SetupHardware();
  void SetupEvaluation();
  void SetupSelection();
  void SetupMutation();
  void SetupDataCollection();

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

  bool ScreenForSolution(prog_org_t & prog_org) {
    eval_util.use_training_set = false;
    begin_program_eval.Trigger(prog_org);
    for (eval_util.current_testID = 0; eval_util.current_testID < TESTING_SET_SIZE; ++eval_util.current_testID) {
      begin_program_test.Trigger(prog_org);
      do_program_test.Trigger(prog_org);
      end_program_test.Trigger(prog_org);
      TestResult result = CalcProgramResultOnTest(prog_org);
      if (!result.pass) {
        end_program_eval.Trigger(prog_org);
        return false;
      }
    }
    end_program_eval.Trigger(prog_org);
    return true;
  }

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
      eval_hardware.Delete();
      inst_lib.Delete();
      prog_world.Delete();
      random.Delete();
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
  SetupProblem();

  // Setup program evaluation.
  std::cout << "==== EXPERIMENT SETUP => evaluation ====" << std::endl;
  SetupEvaluation();

  // Setup program selection.
  std::cout << "==== EXPERIMENT SETUP => selection ====" << std::endl;
  SetupSelection(); 
  
  // Setup program fitness calculations.
  std::cout << "==== EXPERIMENT SETUP => world fitness function (not used by lexicase selection) ====" << std::endl;
  prog_world->SetFitFun([this](prog_org_t & prog_org) {
    double fitness = prog_org.GetPhenotype().total_score;
     if (prog_org.GetPhenotype().num_passes == TRAINING_SET_SIZE) { // Add 'smallness' bonus.
      fitness += ((double)(MAX_PROG_SIZE - prog_org.GetGenome().GetSize()))/(double)MAX_PROG_SIZE;
    }
    return fitness;   
  });

  // Setup data collection
  #ifndef EMSCRIPTEN
  // If we're not compiling to javascript, setup data collection for programs
  std::cout << "==== EXPERIMENT SETUP => data collection ====" << std::endl;
  SetupDataCollection();
  #endif

  // Trigger end of setup.
  std::cout << "==== EXPERIMENT SETUP => triggering end setup signal ====" << std::endl;
  end_setup_sig.Trigger();

  std::cout << "==== EXPERIMENT SETUP => DONE! ====" << std::endl;
  setup = true;
}

/// Run the experiment start->finish [update=0 : update=config.GENERATIONS].
void ProgramSynthesisExperiment::Run() {
  // For each generation, advance 'time' by one step.
  for (update = 0; update <= GENERATIONS; ++update) {
    RunStep();
  }
}

/// Run a single step of the experiment
void ProgramSynthesisExperiment::RunStep() {
  // std::cout << "-- Doing Evaluation --" << std::endl;
  do_evaluation_sig.Trigger();  // (1) Evaluate all members of program population.
  // std::cout << "-- Doing Selection --" << std::endl;
  do_selection_sig.Trigger();   // (2) Select who gets to reproduce!
  // std::cout << "-- Doing Update --" << std::endl;
  do_update_sig.Trigger();      // (3) Run update on relevant worlds (population turnover, etc).
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
  AddDefaultInstructions();

}

void ProgramSynthesisExperiment::SetupProblem() {
  emp_assert(emp::Has(problems, PROBLEM), "Unknown problem!", PROBLEM);
  // Big ol' switch statement to select appropriate problem to setup.
  switch (problems.at(PROBLEM).id) {
    case PROBLEM_ID::NumberIO: { SetupProblem_NumberIO(); break; }
    case PROBLEM_ID::SmallOrLarge: { SetupProblem_SmallOrLarge(); break; }
    case PROBLEM_ID::ForLoopIndex: { SetupProblem_ForLoopIndex(); break; }
    case PROBLEM_ID::CompareStringLengths: { SetupProblem_CompareStringLengths(); break; }
    case PROBLEM_ID::DoubleLetters: { SetupProblem_DoubleLetters(); break; }
    case PROBLEM_ID::CollatzNumbers: { SetupProblem_CollatzNumbers(); break; }
    case PROBLEM_ID::ReplaceSpaceWithNewline: { SetupProblem_ReplaceSpaceWithNewline(); break; }
    case PROBLEM_ID::StringDifferences: { SetupProblem_StringDifferences(); break; }
    case PROBLEM_ID::EvenSquares: { SetupProblem_EvenSquares(); break; }
    case PROBLEM_ID::WallisPi: { SetupProblem_WallisPi(); break; }
    case PROBLEM_ID::StringLengthsBackwards: { SetupProblem_StringLengthsBackwards(); break; }
    case PROBLEM_ID::LastIndexOfZero: { SetupProblem_LastIndexOfZero(); break; }
    case PROBLEM_ID::VectorAverage: { SetupProblem_VectorAverage(); break; }
    case PROBLEM_ID::CountOdds: { SetupProblem_CountOdds(); break; }
    case PROBLEM_ID::MirrorImage: { SetupProblem_MirrorImage(); break; }
    case PROBLEM_ID::SuperAnagrams: { SetupProblem_SuperAnagrams(); break; }
    case PROBLEM_ID::SumOfSquares: { SetupProblem_SumOfSquares(); break; }
    case PROBLEM_ID::VectorsSummed: { SetupProblem_VectorsSummed(); break; }
    case PROBLEM_ID::XWordLines: { SetupProblem_XWordLines(); break; }
    case PROBLEM_ID::PigLatin: { SetupProblem_PigLatin(); break; }
    case PROBLEM_ID::NegativeToZero: { SetupProblem_NegativeToZero(); break; }
    case PROBLEM_ID::ScrabbleScore: { SetupProblem_ScrabbleScore(); break; }
    case PROBLEM_ID::Checksum: { SetupProblem_Checksum(); break; }
    case PROBLEM_ID::Digits: { SetupProblem_Digits(); break; }
    case PROBLEM_ID::Grade: { SetupProblem_Grade(); break; }
    case PROBLEM_ID::Median: { SetupProblem_Median(); break; }
    case PROBLEM_ID::Smallest: { SetupProblem_Smallest(); break; }
    case PROBLEM_ID::Syllables: { SetupProblem_Syllables(); break; }
    default: {
      std::cout << "Unknown problem (" << PROBLEM << "). Exiting." << std::endl;
      exit(-1);
    }
  }
}

void ProgramSynthesisExperiment::SetupEvaluation() {
  std::cout << "Setting up evaluation - every program gets evaluated on all tests in training set." << std::endl;

  // Setup program world on-placement response.
  prog_world->OnPlacement([this](size_t pos) {
    // Reset the program phenotype on placement.
    prog_world->GetOrg(pos).GetPhenotype().Reset(TRAINING_SET_SIZE);
  });

  // What should we happen on evaluation?
  do_evaluation_sig.AddAction([this]() {
    // Evaluate each program on entire training set.
    eval_util.use_training_set = true;
    for (eval_util.current_programID = 0; eval_util.current_programID < PROG_POP_SIZE; ++eval_util.current_programID) {
      emp_assert(prog_world->IsOccupied(eval_util.current_programID));
      prog_org_t & prog_org = prog_world->GetOrg(eval_util.current_programID);
      begin_program_eval.Trigger(prog_org);
      for (eval_util.current_testID = 0; eval_util.current_testID < TRAINING_SET_SIZE; ++eval_util.current_testID) {
        // TestResult result = EvaluateOnTrainingCase(prog_org, eval_util.current_testID);
        begin_program_test.Trigger(prog_org);
        do_program_test.Trigger(prog_org);
        end_program_test.Trigger(prog_org);
        TestResult result = CalcProgramResultOnTest(prog_org);

        // Update program organism's phenotype.
        prog_org_phen_t & prog_phen = prog_org.GetPhenotype();
        prog_phen.RecordScore(eval_util.current_testID, result.score);
        prog_phen.RecordPass(eval_util.current_testID, result.pass);
        prog_phen.RecordSubmission(result.sub);
      }
      end_program_eval.Trigger(prog_org);
    }
  });

  // Post-evaluation, find dominant and screen for solutions.
  do_evaluation_sig.AddAction([this]() {
    double cur_best_score = 0;
    for (eval_util.current_programID = 0; eval_util.current_programID < PROG_POP_SIZE; ++eval_util.current_programID) {
      emp_assert(prog_world->IsOccupied(eval_util.current_programID));
      prog_org_t & prog_org = prog_world->GetOrg(eval_util.current_programID);
      const size_t pass_total = prog_org.GetPhenotype().num_passes;
      const double total_score = prog_org.GetPhenotype().total_score;

      // Is this the highest scoring program this generation?
      if (total_score > cur_best_score || eval_util.current_programID == 0) {
        dominant_prog_id = eval_util.current_programID;
        cur_best_score = total_score;
      }

      // At this point, this program has been evaluated on entire testing set.
      // If it passed all testing set test cases (and it's smaller than any solution
      // we've seen so far), we should check to see if its a solution.
      if (pass_total == TRAINING_SET_SIZE && prog_org.GetGenome().GetSize() < smallest_prog_solution_size) {
        if (ScreenForSolution(prog_org)) { // todo - write screen for solution function
          if (!solution_found) { update_first_solution_found = prog_world->GetUpdate(); }
          solution_found = true;
          smallest_prog_solution_size = prog_org.GetGenome().GetSize();
          // todo -> update solutions file!
        }
      }
    }
  });

}

void ProgramSynthesisExperiment::SetupSelection() {
  std::cout << "Setting up lexicase selection for programs!" << std::endl;
  // 1 function for every test score.
  for (size_t i = 0; i < TRAINING_SET_SIZE; ++i) {
    lexicase_prog_fit_set.push_back([i](prog_org_t & prog_org) {
      emp_assert(i < prog_org.GetPhenotype().test_scores.size(), i, prog_org.GetPhenotype().test_scores.size());
      double score = prog_org.GetPhenotype().test_scores[i];
      return score;
    });
  }
  // Add pressure for small size.
  lexicase_prog_fit_set.push_back([this](prog_org_t & prog_org) {
    if (prog_org.GetPhenotype().num_passes == TRAINING_SET_SIZE) {
      return (double)(MAX_PROG_SIZE - prog_org.GetGenome().GetSize());
    }
    return 0.0;
  });
  // Configure do selection signal.
  do_selection_sig.AddAction([this]() {
    emp::LexicaseSelect_NAIVE(*prog_world,
                              lexicase_prog_fit_set,
                              PROG_POP_SIZE,
                              LEXICASE_MAX_FUNS);
  });
}

void ProgramSynthesisExperiment::SetupMutation() {

  // Configure TagLGP mutator.
  prog_mutator.MAX_PROGRAM_LEN = MAX_PROG_SIZE;
  prog_mutator.MIN_PROGRAM_LEN = MIN_PROG_SIZE;

  prog_mutator.MAX_NUMERIC_ARG = MEM_SIZE-1;

  prog_mutator.PER_BIT_FLIP = PROG_MUT__PER_BIT_FLIP;
  prog_mutator.PER_NUMERIC_ARG_SUB = PROG_MUT__PER_NUMERIC_ARG_SUB;
  prog_mutator.PER_INST_SUB = PROG_MUT__PER_INST_SUB;
  prog_mutator.PER_INST_INS = PROG_MUT__PER_INST_INS;
  prog_mutator.PER_INST_DEL = PROG_MUT__PER_INST_DEL;
  prog_mutator.PER_PROG_SLIP = PROG_MUT__PER_PROG_SLIP;

  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      prog_mutator.USE_NUMERIC_ARGUMENTS = false;
      prog_mutator.USE_TAG_ARGUMENTS = true;
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      prog_mutator.USE_NUMERIC_ARGUMENTS = true;
      prog_mutator.USE_TAG_ARGUMENTS = false;
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      prog_mutator.USE_NUMERIC_ARGUMENTS = true;
      prog_mutator.USE_TAG_ARGUMENTS = true;
      break;
    }
  }
  // If we're not using modules, don't allow module duplications/deletions.
  if (!USE_MODULES) {
    prog_mutator.PER_MOD_DUP = 0.0;
    prog_mutator.PER_MOD_DEL = 0.0;
  }

  // Configure world mutation function.
  prog_world->SetMutFun([this](prog_org_t & prog_org, emp::Random & rnd) {
    return prog_mutator.Mutate(rnd, prog_org.GetGenome());
  });
  
  // Set program world to auto mutate
  end_setup_sig.AddAction([this]() {
    prog_world->SetAutoMutate();
  });

}

void ProgramSynthesisExperiment::SetupDataCollection() {
  std::cout << "Todo!" << std::endl;
}


void ProgramSynthesisExperiment::InitConfigs(const ProgramSynthesisConfig & config) {

  SEED = config.SEED();
  GENERATIONS = config.GENERATIONS();
  PROGRAM_ARGUMENT_MODE = config.PROGRAM_ARGUMENT_MODE();
  PROG_POP_SIZE = config.PROG_POP_SIZE();
  PROBLEM = config.PROBLEM();
  BENCHMARK_DATA_DIR = config.BENCHMARK_DATA_DIR();

  USE_MODULES = config.USE_MODULES();
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

  LEXICASE_MAX_FUNS = config.LEXICASE_MAX_FUNS();

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

void ProgramSynthesisExperiment::AddDefaultInstructions() {
  std::cout << "Adding DEFAULT instructions to instruction set." << std::endl;
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

/// Add instructions (expecting tag arguments) that are shared across all problems.
void ProgramSynthesisExperiment::AddDefaultInstructions_TagArgs() {
  std::cout << "Adding default TAG-BASED ARGUMENT instructions." << std::endl;
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

  if (USE_MODULES) {
    inst_lib->AddInst("Call-Tag", hardware_t::Inst_Call, 1, "Call module using A for tag-based reference");
    inst_lib->AddInst("Routine-Tag", hardware_t::Inst_Routine, 1, "Call module as a routine (don't use call stack)");
    if (!inst_lib->IsInst("ModuleDef")) {
      inst_lib->AddInst("ModuleDef", hardware_t::Inst_Nop, 1, "Define module with tag A", {inst_lib_t::InstProperty::MODULE});
    }
    if (!inst_lib->IsInst("Return")) {
      inst_lib->AddInst("Return", hardware_t::Inst_Return, 0, "Return from current routine/call");
    }
  }
}

/// Add instructions (expecting tag arguments) that are shared across all problems.
void ProgramSynthesisExperiment::AddDefaultInstructions_NumArgs() {
  std::cout << "Adding default NUMERIC ARGUMENT instructions." << std::endl;
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

/// Add numeric terminals
void ProgramSynthesisExperiment::AddNumericTerminals(size_t min, size_t max) {
  std::cout << "Adding NUMERIC TERMINAL instructions to instruction set." << std::endl;
  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      AddNumericTerminals_TagArgs(min, max);
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      AddNumericTerminals_NumArgs(min, max);
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      AddNumericTerminals_NumArgs(min, max);
      AddNumericTerminals_TagArgs(min, max);
      break;
    }
    default: {
      std::cout << "Unrecognized PROGRAM_ARGUMENT_MODE (" << PROGRAM_ARGUMENT_MODE << "). Exiting." << std::endl;
      exit(-1); 
    }
  }
}

void ProgramSynthesisExperiment::AddNumericTerminals_TagArgs(size_t min, size_t max) {
  std::cout << "Adding NUMERIC TERMINAL instructions with TAG-BASED ARGUMENTS to instruction set." << std::endl;
  for (size_t i = min; i <= max; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i) + "-Tag",
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
        if (!hw.IsValidMemPos(posA)) return; // Do nothing
        wmem.Set(posA, (double)i);
      }
    );
  }
}

void ProgramSynthesisExperiment::AddNumericTerminals_TagArgs_NoTypeSearch(size_t min, size_t max) {
  // 
  std::cout << "This isn't implemented!" << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::AddNumericTerminals_NumArgs(size_t min, size_t max) {
  std::cout << "Adding NUMERIC TERMINAL instructions with NUMERIC ARGUMENTS to instruction set." << std::endl;
  for (size_t i = min; i <= max; ++i) {
    inst_lib->AddInst("Set-" + emp::to_string(i) + "-Num",
      [i](hardware_t & hw, const inst_t & inst) {
        hardware_t::CallState & state = hw.GetCurCallState();
        hardware_t::Memory & wmem = state.GetWorkingMem();
        
        emp_assert(inst.arg_nums.size() >= 1);
        size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;
        
        wmem.Set(posA, (double)i);
      }
    );
  }
}

// Problem setups
void ProgramSynthesisExperiment::SetupProblem_NumberIO() {
  std::cout << "Setting up problem: NumberIO." << std::endl;

  using prob_input_t = typename ProblemUtilities_NumberIO::input_t;
  using prob_output_t = typename ProblemUtilities_NumberIO::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  // TODO - CONFIGURE problem utilities!
  prob_utils_NumberIO.MAX_ERROR = emp::Abs(PROB_NUMBER_IO__DOUBLE_MAX) * 2;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_NumberIO.GetTrainingSet().LoadTestCases(training_examples_fpath);
  prob_utils_NumberIO.GetTestingSet().LoadTestCases(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_NumberIO.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_NumberIO.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_NumberIO.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 2);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_NumberIO.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_NumberIO.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input.first);
      wmem.Set(1, input.second);

    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_NumberIO.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_NumberIO.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_NumberIO.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_NumberIO.CalcScoreGradient(correct_output, prob_utils_NumberIO.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);

  // todo - add load and submit instructions (all types)
  
}

void ProgramSynthesisExperiment::SetupProblem_SmallOrLarge() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_ForLoopIndex() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_CompareStringLengths() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_DoubleLetters() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_CollatzNumbers() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_ReplaceSpaceWithNewline() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_StringDifferences() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_EvenSquares() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_WallisPi() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_StringLengthsBackwards() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_LastIndexOfZero() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_VectorAverage() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_CountOdds() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_MirrorImage() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_SuperAnagrams() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_SumOfSquares() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_VectorsSummed() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_XWordLines() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_PigLatin() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_NegativeToZero() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_ScrabbleScore() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_Checksum() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_Digits() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_Grade() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_Median() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_Smallest() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_Syllables() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}



#endif