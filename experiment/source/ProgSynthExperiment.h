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

// Classify instructions: 
// Arg types: TAG_ARGS NUM_ARGS NO_ARGS
// Mem search: MEM_TYPE_SEARCHING MEM_TYPE_NO_SEARCHING MEM_TYPE_AGNOSTIC

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

constexpr size_t PROB_LAST_INDEX_OF_ZERO__MAX_VEC_LEN = 50;

constexpr int PROB_VECTORS_SUMMED__MAX_NUM = 1000;

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
  {"string-lengths-backwards", {PROBLEM_ID::StringLengthsBackwards, "training-examples-string-lengths-backwards.csv", "testing-examples-string-lengths-backwards.csv"}},
  {"last-index-of-zero", {PROBLEM_ID::LastIndexOfZero, "training-examples-last-index-of-zero.csv", "testing-examples-last-index-of-zero.csv"}},
  // {"count-odds", {PROBLEM_ID::CountOdds, "training-examples-count-odds.csv", "testing-examples-count-odds.csv"}},
  {"mirror-image", {PROBLEM_ID::MirrorImage, "training-examples-mirror-image.csv", "testing-examples-mirror-image.csv"}},
  {"vectors-summed", {PROBLEM_ID::VectorsSummed, "training-examples-vectors-summed.csv", "testing-examples-vectors-summed.csv"}},
  // {"sum-of-squares", {PROBLEM_ID::SumOfSquares, "training-examples-sum-of-squares.csv", "testing-examples-sum-of-squares.csv"}},
  // {"vector-average", {PROBLEM_ID::VectorAverage, "training-examples-vector-average.csv", "testing-examples-vector-average.csv"}},
  {"median", {PROBLEM_ID::Median, "training-examples-median.csv", "testing-examples-median.csv"}},
  {"smallest", {PROBLEM_ID::Smallest, "training-examples-smallest.csv", "testing-examples-smallest.csv"}},
  {"grade", {PROBLEM_ID::Grade, "training-examples-grade.csv", "testing-examples-grade.csv"}}
};

class ProgramSynthesisExperiment {
public:

  using hardware_t = typename TagLGP::TagLinearGP_TW<TAG_WIDTH>;
  using inst_lib_t = typename TagLGP::InstLib<hardware_t>;
  using inst_t = typename hardware_t::inst_t;
  using inst_prop_t = typename inst_lib_t::InstProperty;

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
  bool PROGRAM_ARGUMENTS_TYPE_SEARCH;
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

  std::string DATA_DIRECTORY;
  size_t SNAPSHOT_INTERVAL;
  size_t SUMMARY_STATS_INTERVAL;

  // Experiment variables
  size_t TRAINING_SET_SIZE;
  size_t TESTING_SET_SIZE;

  bool USES_VECTOR_INSTRUCTIONS;
  bool USES_STRING_INSTRUCTIONS;

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
  // ProblemUtilities_DoubleLetters prob_utils_DoubleLetters;
  // ProblemUtilities_CollatzNumbers prob_utils_CollatzNumbers;
  // ProblemUtilities_ReplaceSpaceWithNewline prob_utils_ReplaceSpaceWithNewline;
  // ProblemUtilities_StringDifferences prob_utils_StringDifferences;
  // ProblemUtilities_EvenSquares prob_utils_EvenSquares;
  // ProblemUtilities_WallisPi prob_utils_WallisPi;
  ProblemUtilities_StringLengthsBackwards prob_utils_StringLengthsBackwards;
  ProblemUtilities_LastIndexOfZero prob_utils_LastIndexOfZero;
  // ProblemUtilities_VectorAverage prob_utils_VectorAverage;
  // ProblemUtilities_CountOdds prob_utils_CountOdds;
  ProblemUtilities_MirrorImage prob_utils_MirrorImage;
  // ProblemUtilities_SuperAnagrams prob_utils_SuperAnagrams;
  ProblemUtilities_SumOfSquares prob_utils_SumOfSquares;
  ProblemUtilities_VectorsSummed prob_utils_VectorsSummed;
  // ProblemUtilities_XWordLines prob_utils_XWordLines;
  // ProblemUtilities_PigLatin prob_utils_PigLatin;
  // ProblemUtilities_NegativeToZero prob_utils_NegativeToZero;
  // ProblemUtilities_ScrabbleScore prob_utils_ScrabbleScore;
  // ProblemUtilities_Checksum prob_utils_Checksum;
  // ProblemUtilities_Digits prob_utils_Digits;
  ProblemUtilities_Grade prob_utils_Grade;
  ProblemUtilities_Median prob_utils_Median;
  ProblemUtilities_Smallest prob_utils_Smallest;
  // ProblemUtilities_Syllables prob_utils_Syllables;

  emp::Ptr<emp::DataFile> solution_file;
  emp::Ptr<emp::DataFile> prog_arg_info_file;
  
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

    struct TestingSetPhenotype {
      emp::vector<double> test_scores;
      double total_score;
      
      emp::vector<bool> test_passes;
      size_t num_passes;
      size_t num_fails;

      size_t total_submissions;

      void Reset(size_t s=0) {
        // test_results.clear();
        // test_results.resize(s, 0);

        test_scores.clear();
        test_scores.resize(s, 0);
        total_score = 0;
        
        test_passes.clear();
        test_passes.resize(s, false);
        num_passes = 0;
        num_fails = 0;

        total_submissions = 0;
      }

      void RecordScore(size_t id, double val) {
        emp_assert(id < test_scores.size());
        total_score += val;
        test_scores[id] = val;
      }

      void RecordPass(size_t id, bool pass) {
        emp_assert(id < test_passes.size());
        if (pass) ++num_passes;
        else ++num_fails;
        test_passes[id] = pass;
      }

      void RecordSubmission(bool sub) {
        total_submissions += (size_t)sub;
      }
    } testingset_phenotype;

    EvalUtil(size_t pID=0, size_t tID=0) : current_programID(pID), current_testID(tID), use_training_set(true) { ; }
  } eval_util;

  // Experiment signals
  emp::Signal<void(void)> do_evaluation_sig;
  emp::Signal<void(void)> do_selection_sig;
  emp::Signal<void(void)> do_update_sig;  

  emp::Signal<void(void)> do_pop_snapshot_sig;

  emp::Signal<void(void)> end_setup_sig;

  // Program evaluation signals
  emp::Signal<void(prog_org_t &)> begin_program_eval;   ///< Begin evaluating a program on a test case set.
  emp::Signal<void(prog_org_t &)> end_program_eval;     ///< Finish evaluating a program on a test case set.

  emp::Signal<void(prog_org_t &)> begin_program_test;   ///< Begin evaluating a program on an individual test.
  emp::Signal<void(prog_org_t &)> do_program_test;      ///< Do single-test evaluation with given program on given test.
  emp::Signal<void(prog_org_t &)> end_program_test;     ///< Finish evaluating a program on an individual test.

  emp::Signal<void(prog_org_t &)> do_program_advance;   ///< Advance virtual hardware by one time step.

  std::function<TestResult(prog_org_t &)> CalcProgramResultOnTest; ///< Evaluate given program on given training case.

  // Program stats collectors
  struct ProgramStatsCollectors {
    // General
    std::function<size_t(void)> get_id;
    // Fitness evaluation stats
    std::function<double(void)> get_fitness;
    std::function<double(void)> get_fitness_eval__total_score;          
    std::function<size_t(void)> get_fitness_eval__num_passes;           
    std::function<size_t(void)> get_fitness_eval__num_fails;            
    std::function<size_t(void)> get_fitness_eval__num_tests;            
    std::function<std::string(void)> get_fitness_eval__passes_by_test;

    std::function<size_t(void)> get_testingset_eval__num_passes;          // - get_validation_eval__num_passes;
    std::function<size_t(void)> get_testingset_eval__num_tests;           // - get_validation_eval__num_tests
    std::function<std::string(void)> get_testingset_eval__passes_by_test; // - get_validation_eval__passes_by_test

    std::function<size_t(void)> get_program_len;
    std::function<size_t(void)> get_tag_arg_inst_cnt;
    std::function<size_t(void)> get_num_arg_inst_cnt;
    std::function<size_t(void)> get_no_arg_inst_cnt;
    std::function<std::string(void)> get_program;
    
  } program_stats;

  struct ArgTypeDistribution {
    size_t num_arg_instructions;
    size_t tag_arg_instructions;
    size_t no_arg_instructions;
    size_t total_instructions;
  };

  ArgTypeDistribution pop_arg_type_distribution;

  std::function<size_t(void)> get_update;
  
  // Internal functions
  void InitConfigs(const ProgramSynthesisConfig & config);
  
  void InitProgPop_Random();

  void SnapshotPrograms();
  void OutputInstLibInstructionTypes();

  void AddDefaultInstructions();
  void AddDefaultInstructions_TagArgs();
  void AddDefaultInstructions_TagArgs_NoTypeSearch();
  void AddDefaultInstructions_NumArgs();
  void AddDefaultInstructions_NumArgs_WithTypeSearch();

  void AddVectorInstructions_TagArgs();
  void AddVectorInstructions_TagArgs_NoTypeSearch();
  void AddVectorInstructions_NumArgs();
  void AddVectorInstructions_NumArgs_WithTypeSearch();

  void AddStringInstructions_TagArgs();
  void AddStringInstructions_TagArgs_NoTypeSearch();
  void AddStringInstructions_NumArgs();
  void AddStringInstructions_NumArgs_WithTypeSearch();

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

  void DoTestingSetValidation(prog_org_t & prog_org) {
    eval_util.use_training_set = false;
    begin_program_eval.Trigger(prog_org);
    eval_util.testingset_phenotype.Reset(TESTING_SET_SIZE); // Reset eval phenotype
    for (eval_util.current_testID = 0; eval_util.current_testID < TESTING_SET_SIZE; ++eval_util.current_testID) {
      begin_program_test.Trigger(prog_org);
      do_program_test.Trigger(prog_org);
      end_program_test.Trigger(prog_org);
      TestResult result = CalcProgramResultOnTest(prog_org);
      // Update eval phenotype
      eval_util.testingset_phenotype.RecordScore(eval_util.current_testID, result.score);
      eval_util.testingset_phenotype.RecordPass(eval_util.current_testID, result.pass);
      eval_util.testingset_phenotype.RecordSubmission(result.sub);
    }
    end_program_eval.Trigger(prog_org);
  }

  void PrintInstructionSet() {
    std::cout << "=========== Instruction Set ===========" << std::endl;
    inst_lib->Print();
    std::cout << "=======================================" << std::endl;
  }

  ArgTypeDistribution CountArgTypes(const prog_org_gen_t & program) {
    ArgTypeDistribution distribution;
    distribution.num_arg_instructions = 0;
    distribution.tag_arg_instructions = 0;
    distribution.no_arg_instructions = 0;
    distribution.total_instructions = program.GetSize();

    for (size_t i = 0; i < program.GetSize(); ++i) {
      bool has_num_args = inst_lib->HasProperty(program[i].id, inst_prop_t::NUM_ARGS);
      bool has_tag_args = inst_lib->HasProperty(program[i].id, inst_prop_t::TAG_ARGS);
      if (has_num_args) distribution.num_arg_instructions++;
      if (has_tag_args) distribution.tag_arg_instructions++;
      if (!has_num_args && !has_tag_args) distribution.no_arg_instructions++;
    }
    return distribution;
  }

  size_t CountTagArgInstructions(const prog_org_gen_t & program) {
    size_t cnt = 0;
    for (size_t i = 0; i < program.GetSize(); ++i) {
      bool has_tag_args = inst_lib->HasProperty(program[i].id, inst_prop_t::TAG_ARGS);
      if (has_tag_args) cnt++;
    }
    return cnt;
  }

  size_t CountNumArgInstructions(const prog_org_gen_t & program) {
    size_t cnt = 0;
    for (size_t i = 0; i < program.GetSize(); ++i) {
      bool has_num_args = inst_lib->HasProperty(program[i].id, inst_prop_t::NUM_ARGS);
      if (has_num_args) cnt++;
    }
    return cnt;
  }

  size_t CountNoArgInstructions(const prog_org_gen_t & program) {
    size_t cnt = 0;
    for (size_t i = 0; i < program.GetSize(); ++i) {
      bool has_num_args = inst_lib->HasProperty(program[i].id, inst_prop_t::NUM_ARGS);
      bool has_tag_args = inst_lib->HasProperty(program[i].id, inst_prop_t::TAG_ARGS);
      if (!has_num_args && !has_tag_args) cnt++;
    }
    return cnt;
  }


  // ------------ Problem-specific instructions ------------
  // _NO_TYPE_SEARCH
  // _WITH_TYPE_SEARCH
  // -- Number IO --
  // Tag-based args
  void Inst_LoadInt_NumberIO__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadDouble_NumberIO__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_NumberIO__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  // Numeric args
  void Inst_LoadInt_NumberIO__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadDouble_NumberIO__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_NumberIO__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  // with/without type searching
  void Inst_SubmitNum_NumberIO__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_NumberIO__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);
  
  // -- SmallOrLarge --
  // Tag-based args
  void Inst_LoadInt_SmallOrLarge__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  // Numeric args
  void Inst_LoadInt_SmallOrLarge__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  // No args
  void Inst_SubmitSmall_SmallOrLarge(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitLarge_SmallOrLarge(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNone_SmallOrLarge(hardware_t & hw, const inst_t & inst);

  // -- ForLoopIndex --
  // Tag-based args
  void Inst_LoadStart_ForLoopIndex__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadEnd_ForLoopIndex__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadStep_ForLoopIndex__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_ForLoopIndex__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  // Numeric args
  void Inst_LoadStart_ForLoopIndex__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadEnd_ForLoopIndex__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadStep_ForLoopIndex__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_ForLoopIndex__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  // with/without type searching
  void Inst_SubmitNum_ForLoopIndex__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_ForLoopIndex__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);


  // -- Grade --
  // Tag-based args
  void Inst_LoadThreshA_Grade__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadThreshB_Grade__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadThreshC_Grade__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadThreshD_Grade__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadGrade_Grade__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  // Numeric args
  void Inst_LoadThreshA_Grade__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadThreshB_Grade__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadThreshC_Grade__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadThreshD_Grade__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadGrade_Grade__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  // No args
  void Inst_SubmitA_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitB_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitC_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitD_Grade(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitF_Grade(hardware_t & hw, const inst_t & inst);

  // -- Median --
  // Tag-based args
  void Inst_LoadNum1_Median__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum2_Median__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum3_Median__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_Median__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  // Numeric args
  void Inst_LoadNum1_Median__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum2_Median__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum3_Median__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_Median__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  // with/without type searching
  void Inst_SubmitNum_Median__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_Median__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);

  // -- Smallest --
  // Tag-based args
  void Inst_LoadNum1_Smallest__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum2_Smallest__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum3_Smallest__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum4_Smallest__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_Smallest__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  // Numeric args
  void Inst_LoadNum1_Smallest__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum2_Smallest__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum3_Smallest__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadNum4_Smallest__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_Smallest__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  // with/without type searching
  void Inst_SubmitNum_Smallest__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_Smallest__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);

  // -- Compare string lengths --
  void Inst_SubmitTrue_CompareStringLengths(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitFalse_CompareStringLengths(hardware_t & hw, const inst_t & inst);
  // Tag-based args
  void Inst_LoadStr1_CompareStringLengths__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadStr2_CompareStringLengths__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadStr3_CompareStringLengths__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_CompareStringLengths__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  // Numeric args
  void Inst_LoadStr1_CompareStringLengths__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadStr2_CompareStringLengths__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadStr3_CompareStringLengths__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_CompareStringLengths__NUM_ARGS(hardware_t & hw, const inst_t & inst);  
  // with/without type searching
  void Inst_SubmitVal_CompareStringLengths__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_CompareStringLengths__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);

  // -- String length backwards --
  // Tag-based args
  void Inst_LoadStrVec_StringLengthsBackwards__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_StringLengthsBackwards__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  // Numeric args
  void Inst_LoadStrVec_StringLengthsBackwards__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_StringLengthsBackwards__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  // with/without type searching
  void Inst_SubmitVal_StringLengthsBackwards__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_StringLengthsBackwards__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);

  // -- Last index of zero --
  // tag-based args
  void Inst_LoadVec_LastIndexOfZero__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_LastIndexOfZero__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  // numeric args
  void Inst_LoadVec_LastIndexOfZero__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_LastIndexOfZero__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  // with/without type searching
  void Inst_SubmitNum_LastIndexOfZero__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitNum_LastIndexOfZero__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);

  // -- Mirror image --
  // no args
  void Inst_SubmitTrue_MirrorImage(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitFalse_MirrorImage(hardware_t & hw, const inst_t & inst);
  // tag args
  void Inst_LoadVec1_MirrorImage__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadVec2_MirrorImage__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_MirrorImage__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  // num args
  void Inst_LoadVec1_MirrorImage__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadVec2_MirrorImage__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_MirrorImage__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  // with and without type searching
  void Inst_SubmitVal_MirrorImage__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVal_MirrorImage__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);

  // -- Vectors summed --  
  // tag args
  void Inst_LoadVec1_VectorsSummed__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadVec2_VectorsSummed__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVec_VectorsSummed__TAG_ARGS(hardware_t & hw, const inst_t & inst);
  // num args
  void Inst_LoadVec1_VectorsSummed__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_LoadVec2_VectorsSummed__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVec_VectorsSummed__NUM_ARGS(hardware_t & hw, const inst_t & inst);
  // with and without type searching
  void Inst_SubmitVec_VectorsSummed__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);
  void Inst_SubmitVec_VectorsSummed__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst);

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
      solution_file.Delete();
      prog_arg_info_file.Delete();
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
  
  USES_VECTOR_INSTRUCTIONS = false;
  USES_STRING_INSTRUCTIONS = false;

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

    if (update % SNAPSHOT_INTERVAL == 0 || update == GENERATIONS) do_pop_snapshot_sig.Trigger();

    if (update % SUMMARY_STATS_INTERVAL == 0 || update == GENERATIONS) {
      // Update population argument type distribution
      pop_arg_type_distribution.num_arg_instructions = 0;
      pop_arg_type_distribution.tag_arg_instructions = 0;
      pop_arg_type_distribution.no_arg_instructions = 0;
      pop_arg_type_distribution.total_instructions = 0;
      for (size_t i = 0; i < prog_world->GetSize(); ++i) {
        auto prog_distribution = CountArgTypes(prog_world->GetOrg(i).GetGenome());
        pop_arg_type_distribution.num_arg_instructions += prog_distribution.num_arg_instructions;
        pop_arg_type_distribution.tag_arg_instructions += prog_distribution.tag_arg_instructions;
        pop_arg_type_distribution.no_arg_instructions += prog_distribution.no_arg_instructions;
        pop_arg_type_distribution.total_instructions += prog_distribution.total_instructions;
      }
      prog_arg_info_file->Update();
    } 

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
  PrintInstructionSet();

  // Setup program evaluation.
  std::cout << "==== EXPERIMENT SETUP => evaluation ====" << std::endl;
  SetupEvaluation();

  // Setup program selection.
  std::cout << "==== EXPERIMENT SETUP => selection ====" << std::endl;
  SetupSelection(); 

  // Setup program mutation.
  std::cout << "==== EXPERIMENT SETUP => mutation ====" << std::endl;
  SetupMutation();
  
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
  OutputInstLibInstructionTypes();
  #endif

  // Verify instruction set is no nonsense!
  bool okay = true;
  for (size_t i = 0; i < inst_lib->GetSize(); ++i) {
    if (PROGRAM_ARGUMENT_MODE == (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY) {
      if ( (!(inst_lib->HasProperty(i, inst_prop_t::TAG_ARGS) || inst_lib->HasProperty(i, inst_prop_t::NO_ARGS))) || (inst_lib->HasProperty(i, inst_prop_t::NUM_ARGS)) ) { 
        std::cout << "Bad ARG TYPE behavior by: " << inst_lib->GetName(i) << std::endl;
        okay = false;
      }
    } else if (PROGRAM_ARGUMENT_MODE == (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY) {
      if ( (!(inst_lib->HasProperty(i, inst_prop_t::NUM_ARGS) || inst_lib->HasProperty(i, inst_prop_t::NO_ARGS))) || (inst_lib->HasProperty(i, inst_prop_t::TAG_ARGS)) ) { 
        std::cout << "Bad ARG TYPE behavior by: " << inst_lib->GetName(i) << std::endl;
        okay = false;
      }
    } else if (PROGRAM_ARGUMENT_MODE == (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH) {
      if (!(inst_lib->HasProperty(i, inst_prop_t::TAG_ARGS) || inst_lib->HasProperty(i, inst_prop_t::NUM_ARGS) || inst_lib->HasProperty(i, inst_prop_t::NO_ARGS))) { 
        std::cout << "Bad ARG TYPE behavior by: " << inst_lib->GetName(i) << std::endl;
        okay = false;
      }
    } 
    if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
      if (!(inst_lib->HasProperty(i, inst_prop_t::MEM_TYPE_SEARCHING) || inst_lib->HasProperty(i, inst_prop_t::MEM_TYPE_AGNOSTIC))) {
        std::cout << "Bad MEM SEARCH behavior by: " << inst_lib->GetName(i) << std::endl;
        okay = false;
      }
    } else { 
      if (!(inst_lib->HasProperty(i, inst_prop_t::MEM_TYPE_NO_SEARCHING) || inst_lib->HasProperty(i, inst_prop_t::MEM_TYPE_AGNOSTIC))) {
        std::cout << "Bad MEM SEARCH behavior by: " << inst_lib->GetName(i) << std::endl;
        okay = false;
      }
    }
  }
  if (!okay) { std::cout << "At least one instruction definition wasn't okay. Exiting." << std::endl; exit(-1); }

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
          solution_file->Update();
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
  std::cout << "Setting up data collection." << std::endl;
  // Make a directory.
  mkdir(DATA_DIRECTORY.c_str(), ACCESSPERMS);
  if (DATA_DIRECTORY.back() != '/') DATA_DIRECTORY += '/';

  // Setup stats collectors.
  // NOTE - these functions make ALL sorts of assumptions about when they'll be called.
  program_stats.get_id = [this]() { return eval_util.current_programID; };
  program_stats.get_fitness = [this]() { return prog_world->CalcFitnessID(eval_util.current_programID);  };
  program_stats.get_fitness_eval__total_score = [this]() { return prog_world->GetOrg(eval_util.current_programID).GetPhenotype().total_score; };
  program_stats.get_fitness_eval__num_passes = [this]() { return prog_world->GetOrg(eval_util.current_programID).GetPhenotype().num_passes; };
  program_stats.get_fitness_eval__num_fails = [this]() { return prog_world->GetOrg(eval_util.current_programID).GetPhenotype().num_fails; };
  program_stats.get_fitness_eval__num_tests = [this]() { return prog_world->GetOrg(eval_util.current_programID).GetPhenotype().test_passes.size(); };
  program_stats.get_fitness_eval__passes_by_test = [this]() { 
    prog_org_t & prog = prog_world->GetOrg(eval_util.current_programID);
    std::string scores = "\"[";
    for (size_t i = 0; i < prog.GetPhenotype().test_passes.size(); ++i) {
      if (i) scores += ",";
      scores += emp::to_string((size_t)prog.GetPhenotype().test_passes[i]);
    }
    scores += "]\"";
    return scores;
  };
  program_stats.get_testingset_eval__num_passes = [this]() { return eval_util.testingset_phenotype.num_passes; };
  program_stats.get_testingset_eval__num_tests = [this]() { return eval_util.testingset_phenotype.test_passes.size(); };
  program_stats.get_testingset_eval__passes_by_test = [this]() { 
    std::string scores = "\"[";
    for (size_t i = 0; i < eval_util.testingset_phenotype.test_passes.size(); ++i) {
      if (i) scores += ",";
      scores += emp::to_string(eval_util.testingset_phenotype.test_passes[i]);
    }
    scores += "]\"";
    return scores; 
  };
  program_stats.get_program_len = [this]() { return prog_world->GetOrg(eval_util.current_programID).GetGenome().GetSize(); };
  program_stats.get_program = [this]() { 
    std::ostringstream stream;
    prog_world->GetOrg(eval_util.current_programID).GetGenome().PrintCSVEntry(stream);
    return stream.str();
  };

  program_stats.get_tag_arg_inst_cnt = [this]() { return this->CountTagArgInstructions(prog_world->GetOrg(eval_util.current_programID).GetGenome()); };
  program_stats.get_num_arg_inst_cnt = [this]() { return this->CountNumArgInstructions(prog_world->GetOrg(eval_util.current_programID).GetGenome()); };
  program_stats.get_no_arg_inst_cnt = [this]() { return this->CountNoArgInstructions(prog_world->GetOrg(eval_util.current_programID).GetGenome()); };

  get_update = [this]() { return prog_world->GetUpdate(); };

  // Setup solution file.
  solution_file = emp::NewPtr<emp::DataFile>(DATA_DIRECTORY + "/solutions.csv");
  solution_file->AddFun(get_update, "update");
  solution_file->AddFun(program_stats.get_id, "program_id");
  solution_file->AddFun(program_stats.get_tag_arg_inst_cnt, "program_tag_arg_inst_cnt");
  solution_file->AddFun(program_stats.get_num_arg_inst_cnt, "program_num_arg_inst_cnt");
  solution_file->AddFun(program_stats.get_no_arg_inst_cnt, "program_no_arg_inst_cnt");
  solution_file->AddFun(program_stats.get_program_len, "program_len");
  solution_file->AddFun(program_stats.get_program, "program");
  solution_file->PrintHeaderKeys();

  // Setup program population argument distribution file
  prog_arg_info_file = emp::NewPtr<emp::DataFile>(DATA_DIRECTORY + "/inst_arg_types.csv");
  prog_arg_info_file->AddFun(get_update, "update");
  prog_arg_info_file->template AddFun<size_t>([this]() { return pop_arg_type_distribution.num_arg_instructions; }, "total_num_arg_instructions");
  prog_arg_info_file->template AddFun<size_t>([this]() { return pop_arg_type_distribution.tag_arg_instructions; }, "total_tag_arg_instructions");
  prog_arg_info_file->template AddFun<size_t>([this]() { return pop_arg_type_distribution.no_arg_instructions; }, "total_no_arg_instructions");
  prog_arg_info_file->template AddFun<size_t>([this]() { return pop_arg_type_distribution.total_instructions; }, "total_instructions");
  prog_arg_info_file->PrintHeaderKeys();

  do_pop_snapshot_sig.AddAction([this]() {
    SnapshotPrograms();
  });

  prog_world->SetupFitnessFile(DATA_DIRECTORY + "program_fitness.csv").SetTimingRepeat(SUMMARY_STATS_INTERVAL);
  
}


void ProgramSynthesisExperiment::InitConfigs(const ProgramSynthesisConfig & config) {

  SEED = config.SEED();
  GENERATIONS = config.GENERATIONS();
  PROGRAM_ARGUMENT_MODE = config.PROGRAM_ARGUMENT_MODE();
  PROGRAM_ARGUMENTS_TYPE_SEARCH = config.PROGRAM_ARGUMENTS_TYPE_SEARCH();
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

  DATA_DIRECTORY = config.DATA_DIRECTORY();
  SNAPSHOT_INTERVAL = config.SNAPSHOT_INTERVAL();
  SUMMARY_STATS_INTERVAL = config.SUMMARY_STATS_INTERVAL();

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

void ProgramSynthesisExperiment::SnapshotPrograms() {
  std::string snapshot_dir = DATA_DIRECTORY + "pop_" + emp::to_string(prog_world->GetUpdate());
  mkdir(snapshot_dir.c_str(), ACCESSPERMS);

  emp::DataFile file(snapshot_dir + "/program_pop_" + emp::to_string((int)prog_world->GetUpdate()) + ".csv");

  // Add functions to data file.
  file.AddFun(program_stats.get_id, "program_id", "Program ID");
  
  file.AddFun(program_stats.get_fitness, "fitness");
  file.AddFun(program_stats.get_fitness_eval__total_score, "total_score__fitness_eval");
  file.AddFun(program_stats.get_fitness_eval__num_passes, "num_passes__fitness_eval");
  file.AddFun(program_stats.get_fitness_eval__num_fails, "num_fails__fitness_eval");
  file.AddFun(program_stats.get_fitness_eval__num_tests, "num_tests__fitness_eval");
  file.AddFun(program_stats.get_fitness_eval__passes_by_test, "passes_by_test__fitness_eval");

  file.AddFun(program_stats.get_testingset_eval__num_passes, "num_passes__testingset_eval");
  file.AddFun(program_stats.get_testingset_eval__num_tests, "num_tests__testingset_eval");
  file.AddFun(program_stats.get_testingset_eval__passes_by_test, "passes_by_test__testingset_eval");

  file.AddFun(program_stats.get_tag_arg_inst_cnt, "program_tag_arg_inst_cnt");
  file.AddFun(program_stats.get_num_arg_inst_cnt, "program_num_arg_inst_cnt");
  file.AddFun(program_stats.get_no_arg_inst_cnt, "program_no_arg_inst_cnt");
  file.AddFun(program_stats.get_program_len, "program_len");
  file.AddFun(program_stats.get_program, "program");

  file.PrintHeaderKeys();

  // For each program in the population, dump the program and anything we want to know about it.
  eval_util.use_training_set = false;
  for (eval_util.current_programID = 0; eval_util.current_programID < prog_world->GetSize(); ++eval_util.current_programID) {
    if (!prog_world->IsOccupied(eval_util.current_programID)) continue;
    DoTestingSetValidation(prog_world->GetOrg(eval_util.current_programID)); // Do validation for program.
    // Update snapshot file
    file.Update();
  }
}

void ProgramSynthesisExperiment::OutputInstLibInstructionTypes() {
  std::string out_dir = DATA_DIRECTORY;
  mkdir(out_dir.c_str(), ACCESSPERMS);

  emp::DataFile file(out_dir + "/inst_lib_arg_distribution.csv");

  // Count stuff
  size_t total_num_arg_instructions = 0;
  size_t total_tag_arg_instructions = 0;
  size_t total_no_arg_instructions = 0;
  size_t total_mem_searching_instructions = 0;
  size_t total_mem_no_searching_instructions = 0;
  size_t total_mem_agnostic_instruction = 0;

  for (size_t i = 0; i < inst_lib->GetSize(); ++i) {

    if (inst_lib->HasProperty(i, inst_prop_t::NUM_ARGS)) total_num_arg_instructions++;;
    if (inst_lib->HasProperty(i, inst_prop_t::TAG_ARGS)) total_tag_arg_instructions++;
    if (inst_lib->HasProperty(i, inst_prop_t::NO_ARGS)) total_no_arg_instructions++;
    if (inst_lib->HasProperty(i, inst_prop_t::MEM_TYPE_SEARCHING)) total_mem_searching_instructions++;
    if (inst_lib->HasProperty(i, inst_prop_t::MEM_TYPE_NO_SEARCHING)) total_mem_no_searching_instructions++;
    if (inst_lib->HasProperty(i, inst_prop_t::MEM_TYPE_AGNOSTIC)) total_mem_agnostic_instruction++;
    
  }

  // total_instructions
  file.template AddFun<size_t>([this]() { return inst_lib->GetSize(); }, "total_instructions");
  // total_num_arg_instructions
  file.template AddFun<size_t>([this, total_num_arg_instructions]() { return total_num_arg_instructions; }, "total_num_arg_instructions");
  // total_tag_arg_instructions
  file.template AddFun<size_t>([this, total_tag_arg_instructions]() { return total_tag_arg_instructions; }, "total_tag_arg_instructions");
  // total_no_arg_instructions
  file.template AddFun<size_t>([this, total_no_arg_instructions]() { return total_no_arg_instructions; }, "total_no_arg_instructions");
  // total_mem_searching_instructions
  file.template AddFun<size_t>([this, total_mem_searching_instructions]() { return total_mem_searching_instructions; }, "total_mem_searching_instructions");
  // total_mem_no_searching_instructions
  file.template AddFun<size_t>([this, total_mem_no_searching_instructions]() { return total_mem_no_searching_instructions; }, "total_mem_no_searching_instructions");
  // total_mem_agnostic_instruction
  file.template AddFun<size_t>([this, total_mem_agnostic_instruction]() { return total_mem_agnostic_instruction; }, "total_mem_agnostic_instruction");
  
  file.PrintHeaderKeys();
  file.Update();
}

void ProgramSynthesisExperiment::AddDefaultInstructions() {
  std::cout << "Adding DEFAULT instructions to instruction set." << std::endl;
  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) { AddDefaultInstructions_TagArgs(); } 
      else { AddDefaultInstructions_TagArgs_NoTypeSearch(); }
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) { AddDefaultInstructions_NumArgs_WithTypeSearch(); }
      else { AddDefaultInstructions_NumArgs(); }
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        AddDefaultInstructions_TagArgs();
        AddDefaultInstructions_NumArgs_WithTypeSearch();
      } else {
        AddDefaultInstructions_TagArgs_NoTypeSearch();
        AddDefaultInstructions_NumArgs();
      }
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
  inst_lib->AddInst("Add-Tag", hardware_t::Inst_Add, 3, "wmemANY[C] = wmemNUM[A] + wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Sub-Tag", hardware_t::Inst_Sub, 3, "wmemANY[C] = wmemNUM[A] - wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Mult-Tag", hardware_t::Inst_Mult, 3, "wmemANY[C] = wmemNUM[A] * wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Div-Tag", hardware_t::Inst_Div, 3, "if (wmemNUM[B] != 0) wmemANY[C] = wmemNUM[A] / wmemNUM[B]; else NOP", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Mod-Tag", hardware_t::Inst_Mod, 3, "if (wmemNUM[B] != 0) wmemANY[C] = int(wmemNUM[A]) % int(wmemNUM[B]); else NOP", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumEqu-Tag", hardware_t::Inst_TestNumEqu, 3, "wmemANY[C] = wmemNUM[A] == wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumNEqu-Tag", hardware_t::Inst_TestNumNEqu, 3, "wmemANY[C] = wmemNUM[A] != wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumLess-Tag", hardware_t::Inst_TestNumLess, 3, "wmemANY[C] = wmemNUM[A] < wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumLessTEqu-Tag", hardware_t::Inst_TestNumLessTEqu, 3, "wmemANY[C] = wmemNUM[A] <= wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumGreater-Tag", hardware_t::Inst_TestNumGreater, 3, "wmemANY[C] = wmemNUM[A] > wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumGreaterTEqu-Tag", hardware_t::Inst_TestNumGreaterTEqu, 3, "wmemANY[C] = wmemNUM[A] >= wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Floor-Tag", hardware_t::Inst_Floor, 1, "wmemNUM[A] = floor(wmemNUM[A])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Not-Tag", hardware_t::Inst_Not, 1, "wmemNUM[A] = !wmemNUM[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING}); 
  inst_lib->AddInst("Inc-Tag", hardware_t::Inst_Inc, 1, "wmemNUM[A] = wmemNUM[A] + 1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Dec-Tag", hardware_t::Inst_Dec, 1, "wmemNUM[A] = wmemNUM[A] - 1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});

  // - Memory-related instructions -
  inst_lib->AddInst("CopyMem-Tag", hardware_t::Inst_CopyMem, 2, "wmemANY[B] = wmemANY[A] // Copy mem[A] to mem[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SwapMem-Tag", hardware_t::Inst_SwapMem, 2, "swap(wmemANY[A], wmemANY[B])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

  // - Non-module flow control instructions -
  inst_lib->AddInst("If-Tag", hardware_t::Inst_If, 1, "Execute next flow if(wmemANY[A]) // To be true, mem loc must be non-zero number", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("IfNot-Tag", hardware_t::Inst_IfNot, 1, "Execute next flow if(!wmemANY[A])", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("While-Tag", hardware_t::Inst_While, 1, "While loop over wmemANY[A]", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("Countdown-Tag", hardware_t::Inst_Countdown, 1, "Countdown loop with wmemANY as index.", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  
  // The below instructions take no arguments, check if instruction in library first.
  if (!inst_lib->IsInst("Close")) {
    inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "Close flow", {inst_lib_t::InstProperty::END_FLOW, inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
  if (!inst_lib->IsInst("Break")) {
    inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "Break current flow", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
  if (!inst_lib->IsInst("Return")) {
    inst_lib->AddInst("Return", hardware_t::Inst_Return, 0, "Return from current routine/call", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }

  // Add multi-type instructions (instructions that are only necessary if using multiple types)
  if (USES_VECTOR_INSTRUCTIONS) {
    AddVectorInstructions_TagArgs();
  }
  if (USES_STRING_INSTRUCTIONS) {
    AddStringInstructions_TagArgs();
  }
  if (USES_STRING_INSTRUCTIONS || USES_VECTOR_INSTRUCTIONS) {
    inst_lib->AddInst("IsNum-Tag", hardware_t::Inst_IsNum, 2, "wmemANY[B] = IsNum(wmemANY[A])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
    inst_lib->AddInst("TestMemEqu-Tag", hardware_t::Inst_TestMemEqu, 3, "wmemANY[C] = wmemANY[A] == wmemANY[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
    inst_lib->AddInst("TestMemNEqu-Tag", hardware_t::Inst_TestMemNEqu, 3, "wmemANY[C] = wmemANY[A] != wmemANY[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
}

/// Add instructions (expecting tag arguments)
void ProgramSynthesisExperiment::AddDefaultInstructions_TagArgs_NoTypeSearch() {
  std::cout << "Adding default TAG-BASED ARGUMENT instructions WITHOUT TYPE SEARCHING." << std::endl;
  // - Numeric-related instructions -
  inst_lib->AddInst("Add-Tag", hardware_t::Inst_Add__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] + wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Sub-Tag", hardware_t::Inst_Sub__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] - wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Mult-Tag", hardware_t::Inst_Mult__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] * wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Div-Tag", hardware_t::Inst_Div__TAG_ARGS_NO_TYPE_SEARCH, 3, "if (wmemNUM[B] != 0) wmemANY[C] = wmemNUM[A] / wmemNUM[B]; else NOP", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Mod-Tag", hardware_t::Inst_Mod__TAG_ARGS_NO_TYPE_SEARCH, 3, "if (wmemNUM[B] != 0) wmemANY[C] = int(wmemNUM[A]) % int(wmemNUM[B]); else NOP", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumEqu-Tag", hardware_t::Inst_TestNumEqu__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] == wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumNEqu-Tag", hardware_t::Inst_TestNumNEqu__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] != wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumLess-Tag", hardware_t::Inst_TestNumLess__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] < wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumLessTEqu-Tag", hardware_t::Inst_TestNumLessTEqu__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] <= wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumGreater-Tag", hardware_t::Inst_TestNumGreater__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] > wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumGreaterTEqu-Tag", hardware_t::Inst_TestNumGreaterTEqu__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] >= wmemNUM[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Floor-Tag", hardware_t::Inst_Floor__TAG_ARGS_NO_TYPE_SEARCH, 1, "wmemNUM[A] = floor(wmemNUM[A])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Not-Tag", hardware_t::Inst_Not__TAG_ARGS_NO_TYPE_SEARCH, 1, "wmemNUM[A] = !wmemNUM[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING}); 
  inst_lib->AddInst("Inc-Tag", hardware_t::Inst_Inc__TAG_ARGS_NO_TYPE_SEARCH, 1, "wmemNUM[A] = wmemNUM[A] + 1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Dec-Tag", hardware_t::Inst_Dec__TAG_ARGS_NO_TYPE_SEARCH, 1, "wmemNUM[A] = wmemNUM[A] - 1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});

  // - Memory-related instructions - (these don't check type)
  inst_lib->AddInst("CopyMem-Tag", hardware_t::Inst_CopyMem, 2, "wmemANY[B] = wmemANY[A] // Copy mem[A] to mem[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SwapMem-Tag", hardware_t::Inst_SwapMem, 2, "swap(wmemANY[A], wmemANY[B])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

  // - Non-module flow control instructions - (these don't check type)
  inst_lib->AddInst("If-Tag", hardware_t::Inst_If, 1, "Execute next flow if(wmemANY[A]) // To be true, mem loc must be non-zero number", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("IfNot-Tag", hardware_t::Inst_IfNot, 1, "Execute next flow if(!wmemANY[A])", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("While-Tag", hardware_t::Inst_While, 1, "While loop over wmemANY[A]", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("Countdown-Tag", hardware_t::Inst_Countdown, 1, "Countdown loop with wmemANY as index.", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  
  // The below instructions take no arguments, check if instruction in library first.
  if (!inst_lib->IsInst("Close")) {
    inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "Close flow", {inst_lib_t::InstProperty::END_FLOW, inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
  if (!inst_lib->IsInst("Break")) {
    inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "Break current flow", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
  if (!inst_lib->IsInst("Return")) {
    inst_lib->AddInst("Return", hardware_t::Inst_Return, 0, "Return from current routine/call", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }

  // Add multi-type instructions (instructions that are only necessary if using multiple types)
  if (USES_VECTOR_INSTRUCTIONS) {
    AddVectorInstructions_TagArgs_NoTypeSearch();
  }
  if (USES_STRING_INSTRUCTIONS) {
    AddStringInstructions_TagArgs_NoTypeSearch();
  }
  if (USES_STRING_INSTRUCTIONS || USES_VECTOR_INSTRUCTIONS) {
    inst_lib->AddInst("IsNum-Tag", hardware_t::Inst_IsNum, 2, "wmemANY[B] = IsNum(wmemANY[A])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
    inst_lib->AddInst("TestMemEqu-Tag", hardware_t::Inst_TestMemEqu, 3, "wmemANY[C] = wmemANY[A] == wmemANY[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
    inst_lib->AddInst("TestMemNEqu-Tag", hardware_t::Inst_TestMemNEqu, 3, "wmemANY[C] = wmemANY[A] != wmemANY[B]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
}

/// Add instructions (expecting numeric arguments) that are shared across all problems.
void ProgramSynthesisExperiment::AddDefaultInstructions_NumArgs() {
  std::cout << "Adding default NUMERIC ARGUMENT instructions." << std::endl;
  // - Numeric-related instructions -
  inst_lib->AddInst("Add-Num", hardware_t::Inst_Add__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] + wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Sub-Num", hardware_t::Inst_Sub__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] - wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Mult-Num", hardware_t::Inst_Mult__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] * wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Div-Num", hardware_t::Inst_Div__NUM_ARGS, 3, "if (wmemNUM[B] != 0) wmemANY[C] = wmemNUM[A] / wmemNUM[B]; else NOP", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Mod-Num", hardware_t::Inst_Mod__NUM_ARGS, 3, "if (wmemNUM[B] != 0) wmemANY[C] = int(wmemNUM[A]) % int(wmemNUM[B]); else NOP", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumEqu-Num", hardware_t::Inst_TestNumEqu__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] == wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumNEqu-Num", hardware_t::Inst_TestNumNEqu__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] != wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumLess-Num", hardware_t::Inst_TestNumLess__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] < wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumLessTEqu-Num", hardware_t::Inst_TestNumLessTEqu__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] <= wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumGreater-Num", hardware_t::Inst_TestNumGreater__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] > wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("TestNumGreaterTEqu-Num", hardware_t::Inst_TestNumGreaterTEqu__NUM_ARGS, 3, "wmemANY[C] = wmemNUM[A] >= wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Floor-Num", hardware_t::Inst_Floor__NUM_ARGS, 1, "wmemNUM[A] = floor(wmemNUM[A])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Not-Num", hardware_t::Inst_Not__NUM_ARGS, 1, "wmemNUM[A] = !wmemNUM[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING}); 
  inst_lib->AddInst("Inc-Num", hardware_t::Inst_Inc__NUM_ARGS, 1, "wmemNUM[A] = wmemNUM[A] + 1", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Dec-Num", hardware_t::Inst_Dec__NUM_ARGS, 1, "wmemNUM[A] = wmemNUM[A] - 1", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});

  // - Memory-related instructions -
  inst_lib->AddInst("CopyMem-Num", hardware_t::Inst_CopyMem__NUM_ARGS, 2, "wmemANY[B] = wmemANY[A] // Copy mem[A] to mem[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SwapMem-Num", hardware_t::Inst_SwapMem__NUM_ARGS, 2, "swap(wmemANY[A], wmemANY[B])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

  // - Non-module flow control instructions -
  inst_lib->AddInst("If-Num", hardware_t::Inst_If__NUM_ARGS, 1, "Execute next flow if(wmemANY[A]) // To be true, mem loc must be non-zero number", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("IfNot-Num", hardware_t::Inst_IfNot__NUM_ARGS, 1, "Execute next flow if(!wmemANY[A])", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("While-Num", hardware_t::Inst_While__NUM_ARGS, 1, "While loop over wmemANY[A]", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("Countdown-Num", hardware_t::Inst_Countdown__NUM_ARGS, 1, "Countdown loop with wmemANY as index.", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  
  // The below instructions take no arguments, check if instruction in library first.
  if (!inst_lib->IsInst("Close")) {
    inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "Close flow", {inst_lib_t::InstProperty::END_FLOW, inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
  if (!inst_lib->IsInst("Break")) {
    inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "Break current flow", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
  if (!inst_lib->IsInst("Return")) {
    inst_lib->AddInst("Return", hardware_t::Inst_Return, 0, "Return from current routine/call", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }

  // Add multi-type instructions (instructions that are only necessary if using multiple types)
  if (USES_VECTOR_INSTRUCTIONS) {
    AddVectorInstructions_NumArgs();
  }
  if (USES_STRING_INSTRUCTIONS) {
    AddStringInstructions_NumArgs();
  }
  if (USES_STRING_INSTRUCTIONS || USES_VECTOR_INSTRUCTIONS) {
    inst_lib->AddInst("IsNum-Num", hardware_t::Inst_IsNum__NUM_ARGS, 2, "wmemANY[B] = IsNum(wmemANY[A])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
    inst_lib->AddInst("TestMemEqu-Num", hardware_t::Inst_TestMemEqu__NUM_ARGS, 3, "wmemANY[C] = wmemANY[A] == wmemANY[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
    inst_lib->AddInst("TestMemNEqu-Num", hardware_t::Inst_TestMemNEqu__NUM_ARGS, 3, "wmemANY[C] = wmemANY[A] != wmemANY[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
}

/// Add instructions
void ProgramSynthesisExperiment::AddDefaultInstructions_NumArgs_WithTypeSearch() {
  std::cout << "Adding default NUMERIC ARGUMENT instructions." << std::endl;
  // - Numeric-related instructions -
  inst_lib->AddInst("Add-Num", hardware_t::Inst_Add__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] + wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Sub-Num", hardware_t::Inst_Sub__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] - wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Mult-Num", hardware_t::Inst_Mult__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] * wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Div-Num", hardware_t::Inst_Div__NUM_ARGS_WITH_TYPE_SEARCH, 3, "if (wmemNUM[B] != 0) wmemANY[C] = wmemNUM[A] / wmemNUM[B]; else NOP", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Mod-Num", hardware_t::Inst_Mod__NUM_ARGS_WITH_TYPE_SEARCH, 3, "if (wmemNUM[B] != 0) wmemANY[C] = int(wmemNUM[A]) % int(wmemNUM[B]); else NOP", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumEqu-Num", hardware_t::Inst_TestNumEqu__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] == wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumNEqu-Num", hardware_t::Inst_TestNumNEqu__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] != wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumLess-Num", hardware_t::Inst_TestNumLess__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] < wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumLessTEqu-Num", hardware_t::Inst_TestNumLessTEqu__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] <= wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumGreater-Num", hardware_t::Inst_TestNumGreater__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] > wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("TestNumGreaterTEqu-Num", hardware_t::Inst_TestNumGreaterTEqu__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C] = wmemNUM[A] >= wmemNUM[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Floor-Num", hardware_t::Inst_Floor__NUM_ARGS_WITH_TYPE_SEARCH, 1, "wmemNUM[A] = floor(wmemNUM[A])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Not-Num", hardware_t::Inst_Not__NUM_ARGS_WITH_TYPE_SEARCH, 1, "wmemNUM[A] = !wmemNUM[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Inc-Num", hardware_t::Inst_Inc__NUM_ARGS_WITH_TYPE_SEARCH, 1, "wmemNUM[A] = wmemNUM[A] + 1", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Dec-Num", hardware_t::Inst_Dec__NUM_ARGS_WITH_TYPE_SEARCH, 1, "wmemNUM[A] = wmemNUM[A] - 1", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});

  // - Memory-related instructions - (type-agnostic)
  inst_lib->AddInst("CopyMem-Num", hardware_t::Inst_CopyMem__NUM_ARGS, 2, "wmemANY[B] = wmemANY[A] // Copy mem[A] to mem[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SwapMem-Num", hardware_t::Inst_SwapMem__NUM_ARGS, 2, "swap(wmemANY[A], wmemANY[B])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

  // - Non-module flow control instructions - (type-agnostic)
  inst_lib->AddInst("If-Num", hardware_t::Inst_If__NUM_ARGS, 1, "Execute next flow if(wmemANY[A]) // To be true, mem loc must be non-zero number", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("IfNot-Num", hardware_t::Inst_IfNot__NUM_ARGS, 1, "Execute next flow if(!wmemANY[A])", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("While-Num", hardware_t::Inst_While__NUM_ARGS, 1, "While loop over wmemANY[A]", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("Countdown-Num", hardware_t::Inst_Countdown__NUM_ARGS, 1, "Countdown loop with wmemANY as index.", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  
  // The below instructions take no arguments, check if instruction in library first.
  if (!inst_lib->IsInst("Close")) {
    inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "Close flow", {inst_lib_t::InstProperty::END_FLOW, inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
  if (!inst_lib->IsInst("Break")) {
    inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "Break current flow", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
  if (!inst_lib->IsInst("Return")) {
    inst_lib->AddInst("Return", hardware_t::Inst_Return, 0, "Return from current routine/call", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }

  // Add multi-type instructions (instructions that are only necessary if using multiple types)
  if (USES_VECTOR_INSTRUCTIONS) {
    AddVectorInstructions_NumArgs_WithTypeSearch();
  }
  if (USES_STRING_INSTRUCTIONS) {
    AddStringInstructions_NumArgs_WithTypeSearch();
  }
  if (USES_STRING_INSTRUCTIONS || USES_VECTOR_INSTRUCTIONS) {
    inst_lib->AddInst("IsNum-Num", hardware_t::Inst_IsNum__NUM_ARGS, 2, "wmemANY[B] = IsNum(wmemANY[A])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
    inst_lib->AddInst("TestMemEqu-Num", hardware_t::Inst_TestMemEqu__NUM_ARGS, 3, "wmemANY[C] = wmemANY[A] == wmemANY[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
    inst_lib->AddInst("TestMemNEqu-Num", hardware_t::Inst_TestMemNEqu__NUM_ARGS, 3, "wmemANY[C] = wmemANY[A] != wmemANY[B]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
}

/// Add vector instructions (expecting tag arguments)
void ProgramSynthesisExperiment::AddVectorInstructions_TagArgs() {
  // Memory type agnostic
  inst_lib->AddInst("IsVec-Tag", hardware_t::Inst_IsVec, 2, "wmemANY[B] = IsVec(wmemANY[A])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("MakeVector-Tag", hardware_t::Inst_MakeVector, 3, "wmemANY[C]=Vector([wmemANY[min(A,B),max(A,B)])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  // Memory type matters
  inst_lib->AddInst("VecGet-Tag", hardware_t::Inst_VecGet, 3, "wmemANY[C]=wmemVEC[A][wmemNUM[B]]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecSet-Tag", hardware_t::Inst_VecSet, 3, "wmemVEC[A][wmemNUM[B]]=wmemNUM/STR[C]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecLen-Tag", hardware_t::Inst_VecLen, 2, "wmemANY[B]=wmemVEC[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecAppend-Tag", hardware_t::Inst_VecAppend, 2, "wmemVEC[A].Append(wmemNUM/STR[B])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecPop-Tag", hardware_t::Inst_VecPop, 2, "wmemANY[B]=wmemVEC[A].pop()", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecRemove-Tag", hardware_t::Inst_VecRemove, 2, "wmemVEC[A].Remove(wmemNUM[B])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecReplaceAll-Tag", hardware_t::Inst_VecReplaceAll, 3, "Replace all values (wmemNUM/STR[B]) in wmemVEC[A] with wmemNUM/STR[C]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecIndexOf-Tag", hardware_t::Inst_VecIndexOf, 3, "wmemANY[C] = index of wmemNUM/STR[B] in wmemVEC[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecOccurrencesOf-Tag", hardware_t::Inst_VecOccurrencesOf, 3, "wmemANY[C]= occurrances of wmemNUM/STR[B] in wmemVEC[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecReverse-Tag", hardware_t::Inst_VecReverse, 1, "wmemVEC[A] = Reverse(wmemVEC[A])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecSwapIfLess-Tag", hardware_t::Inst_VecSwapIfLess, 3, "Swap two indices in wmemVEC[A] if vec[wmemNUM[A]] < vec[wmemNUM[B]].", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecGetFront-Tag", hardware_t::Inst_VecGetFront, 2, "wmemANY[B] = front of wmemVEC[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecGetBack-Tag", hardware_t::Inst_VecGetBack, 2, "wmemANY[B] = back of wmemVEC[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Foreach-Tag", hardware_t::Inst_Foreach, 2, "For each thing in wmemVEC[B]", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
}

void ProgramSynthesisExperiment::AddVectorInstructions_TagArgs_NoTypeSearch() {
  // Memory type agnostic
  inst_lib->AddInst("IsVec-Tag", hardware_t::Inst_IsVec, 2, "wmemANY[B] = IsVec(wmemANY[A])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("MakeVector-Tag", hardware_t::Inst_MakeVector, 3, "wmemANY[C]=Vector([wmemANY[min(A,B),max(A,B)])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  // Memory type matters
  inst_lib->AddInst("VecGet-Tag", hardware_t::Inst_VecGet__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C]=wmemVEC[A][wmemNUM[B]]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecSet-Tag", hardware_t::Inst_VecSet__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemVEC[A][wmemNUM[B]]=wmemNUM/STR[C]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecLen-Tag", hardware_t::Inst_VecLen__TAG_ARGS_NO_TYPE_SEARCH, 2, "wmemANY[B]=wmemVEC[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecAppend-Tag", hardware_t::Inst_VecAppend__TAG_ARGS_NO_TYPE_SEARCH, 2, "wmemVEC[A].Append(wmemNUM/STR[B])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecPop-Tag", hardware_t::Inst_VecPop__TAG_ARGS_NO_TYPE_SEARCH, 2, "wmemANY[B]=wmemVEC[A].pop()", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecRemove-Tag", hardware_t::Inst_VecRemove__TAG_ARGS_NO_TYPE_SEARCH, 2, "wmemVEC[A].Remove(wmemNUM[B])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecReplaceAll-Tag", hardware_t::Inst_VecReplaceAll__TAG_ARGS_NO_TYPE_SEARCH, 3, "Replace all values (wmemNUM/STR[B]) in wmemVEC[A] with wmemNUM/STR[C]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecIndexOf-Tag", hardware_t::Inst_VecIndexOf__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C] = index of wmemNUM/STR[B] in wmemVEC[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecOccurrencesOf-Tag", hardware_t::Inst_VecOccurrencesOf__TAG_ARGS_NO_TYPE_SEARCH, 3, "wmemANY[C]= occurrances of wmemNUM/STR[B] in wmemVEC[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecReverse-Tag", hardware_t::Inst_VecReverse__TAG_ARGS_NO_TYPE_SEARCH, 1, "wmemVEC[A] = Reverse(wmemVEC[A])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecSwapIfLess-Tag", hardware_t::Inst_VecSwapIfLess__TAG_ARGS_NO_TYPE_SEARCH, 3, "Swap two indices in wmemVEC[A] if vec[wmemNUM[A]] < vec[wmemNUM[B]].", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecGetFront-Tag", hardware_t::Inst_VecGetFront__TAG_ARGS_NO_TYPE_SEARCH, 2, "wmemANY[B] = front of wmemVEC[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecGetBack-Tag", hardware_t::Inst_VecGetBack__TAG_ARGS_NO_TYPE_SEARCH, 2, "wmemANY[B] = back of wmemVEC[A]", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Foreach-Tag", hardware_t::Inst_Foreach__TAG_ARGS_NO_TYPE_SEARCH, 2, "For each thing in wmemVEC[B]", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});

}

/// Add vector-specific instructions (expecting numeric arguments)
void ProgramSynthesisExperiment::AddVectorInstructions_NumArgs() {
  inst_lib->AddInst("IsVec-Num", hardware_t::Inst_IsVec__NUM_ARGS, 2, "wmemANY[B] = IsVec(wmemANY[A])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("MakeVector-Num", hardware_t::Inst_MakeVector__NUM_ARGS, 3, "wmemANY[C]=Vector([wmemANY[min(A,B),max(A,B)])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

  inst_lib->AddInst("VecGet-Num", hardware_t::Inst_VecGet__NUM_ARGS, 3, "wmemANY[C]=wmemVEC[A][wmemNUM[B]]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecSet-Num", hardware_t::Inst_VecSet__NUM_ARGS, 3, "wmemVEC[A][wmemNUM[B]]=wmemNUM/STR[C]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecLen-Num", hardware_t::Inst_VecLen__NUM_ARGS, 2, "wmemANY[B]=wmemVEC[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecAppend-Num", hardware_t::Inst_VecAppend__NUM_ARGS, 2, "wmemVEC[A].Append(wmemNUM/STR[B])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecPop-Num", hardware_t::Inst_VecPop__NUM_ARGS, 2, "wmemANY[B]=wmemVEC[A].pop()", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecRemove-Num", hardware_t::Inst_VecRemove__NUM_ARGS, 2, "wmemVEC[A].Remove(wmemNUM[B])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecReplaceAll-Num", hardware_t::Inst_VecReplaceAll__NUM_ARGS, 3, "Replace all values (wmemNUM/STR[B]) in wmemVEC[A] with wmemNUM/STR[C]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecIndexOf-Num", hardware_t::Inst_VecIndexOf__NUM_ARGS, 3, "wmemANY[C] = index of wmemNUM/STR[B] in wmemVEC[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecOccurrencesOf-Num", hardware_t::Inst_VecOccurrencesOf__NUM_ARGS, 3, "wmemANY[C]= occurrances of wmemNUM/STR[B] in wmemVEC[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecReverse-Num", hardware_t::Inst_VecReverse__NUM_ARGS, 1, "wmemVEC[A] = Reverse(wmemVEC[A])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecSwapIfLess-Num", hardware_t::Inst_VecSwapIfLess__NUM_ARGS, 3, "Swap two indices in wmemVEC[A] if vec[wmemNUM[A]] < vec[wmemNUM[B]].", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecGetFront-Num", hardware_t::Inst_VecGetFront__NUM_ARGS, 2, "wmemANY[B] = front of wmemVEC[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("VecGetBack-Num", hardware_t::Inst_VecGetBack__NUM_ARGS, 2, "wmemANY[B] = back of wmemVEC[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("Foreach-Num", hardware_t::Inst_Foreach__NUM_ARGS, 2, "For each thing in wmemVEC[B]", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
}

void ProgramSynthesisExperiment::AddVectorInstructions_NumArgs_WithTypeSearch() {
  inst_lib->AddInst("IsVec-Num", hardware_t::Inst_IsVec__NUM_ARGS, 2, "wmemANY[B] = IsVec(wmemANY[A])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("MakeVector-Num", hardware_t::Inst_MakeVector__NUM_ARGS, 3, "wmemANY[C]=Vector([wmemANY[min(A,B),max(A,B)])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 

  inst_lib->AddInst("VecGet-Num", hardware_t::Inst_VecGet__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C]=wmemVEC[A][wmemNUM[B]]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecSet-Num", hardware_t::Inst_VecSet__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemVEC[A][wmemNUM[B]]=wmemNUM/STR[C]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecLen-Num", hardware_t::Inst_VecLen__NUM_ARGS_WITH_TYPE_SEARCH, 2, "wmemANY[B]=wmemVEC[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecAppend-Num", hardware_t::Inst_VecAppend__NUM_ARGS_WITH_TYPE_SEARCH, 2, "wmemVEC[A].Append(wmemNUM/STR[B])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecPop-Num", hardware_t::Inst_VecPop__NUM_ARGS_WITH_TYPE_SEARCH, 2, "wmemANY[B]=wmemVEC[A].pop()", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecRemove-Num", hardware_t::Inst_VecRemove__NUM_ARGS_WITH_TYPE_SEARCH, 2, "wmemVEC[A].Remove(wmemNUM[B])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecReplaceAll-Num", hardware_t::Inst_VecReplaceAll__NUM_ARGS_WITH_TYPE_SEARCH, 3, "Replace all values (wmemNUM/STR[B]) in wmemVEC[A] with wmemNUM/STR[C]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecIndexOf-Num", hardware_t::Inst_VecIndexOf__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C] = index of wmemNUM/STR[B] in wmemVEC[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecOccurrencesOf-Num", hardware_t::Inst_VecOccurrencesOf__NUM_ARGS_WITH_TYPE_SEARCH, 3, "wmemANY[C]= occurrances of wmemNUM/STR[B] in wmemVEC[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecReverse-Num", hardware_t::Inst_VecReverse__NUM_ARGS_WITH_TYPE_SEARCH, 1, "wmemVEC[A] = Reverse(wmemVEC[A])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecSwapIfLess-Num", hardware_t::Inst_VecSwapIfLess__NUM_ARGS_WITH_TYPE_SEARCH, 3, "Swap two indices in wmemVEC[A] if vec[wmemNUM[A]] < vec[wmemNUM[B]].", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecGetFront-Num", hardware_t::Inst_VecGetFront__NUM_ARGS_WITH_TYPE_SEARCH, 2, "wmemANY[B] = front of wmemVEC[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("VecGetBack-Num", hardware_t::Inst_VecGetBack__NUM_ARGS_WITH_TYPE_SEARCH, 2, "wmemANY[B] = back of wmemVEC[A]", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("Foreach-Num", hardware_t::Inst_Foreach__NUM_ARGS_WITH_TYPE_SEARCH, 2, "For each thing in wmemVEC[B]", {inst_lib_t::InstProperty::BEGIN_FLOW, inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
}

void ProgramSynthesisExperiment::AddStringInstructions_TagArgs() {
  inst_lib->AddInst("IsStr-Tag", hardware_t::Inst_IsStr, 2, "wmemANY[B] = IsStr(wmemANY[A])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

  inst_lib->AddInst("StrLength-Tag", hardware_t::Inst_StrLength, 2, "StrLength", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("StrConcat-Tag", hardware_t::Inst_StrConcat, 3, "StrConcat", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
}

void ProgramSynthesisExperiment::AddStringInstructions_TagArgs_NoTypeSearch() {
  inst_lib->AddInst("IsStr-Tag", hardware_t::Inst_IsStr, 2, "wmemANY[B] = IsStr(wmemANY[A])", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

  inst_lib->AddInst("StrLength-Tag", hardware_t::Inst_StrLength__TAG_ARGS_NO_TYPE_SEARCH, 2, "StrLength", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("StrConcat-Tag", hardware_t::Inst_StrConcat__TAG_ARGS_NO_TYPE_SEARCH, 3, "StrConcat", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
}

void ProgramSynthesisExperiment::AddStringInstructions_NumArgs() {
  inst_lib->AddInst("IsStr-Num", hardware_t::Inst_IsStr__NUM_ARGS, 2, "wmemANY[B] = IsStr(wmemANY[A])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

  inst_lib->AddInst("StrLength-Num", hardware_t::Inst_StrLength__NUM_ARGS, 2, "StrLength", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
  inst_lib->AddInst("StrConcat-Num", hardware_t::Inst_StrConcat__NUM_ARGS, 3, "StrConcat", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
}

void ProgramSynthesisExperiment::AddStringInstructions_NumArgs_WithTypeSearch() {
  inst_lib->AddInst("IsStr-Num", hardware_t::Inst_IsStr__NUM_ARGS, 2, "wmemANY[B] = IsStr(wmemANY[A])", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

  inst_lib->AddInst("StrLength-Num", hardware_t::Inst_StrLength__NUM_ARGS_WITH_TYPE_SEARCH, 2, "StrLength", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
  inst_lib->AddInst("StrConcat-Num", hardware_t::Inst_StrConcat__NUM_ARGS_WITH_TYPE_SEARCH, 3, "StrConcat", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
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
      }, 1, "Terminal value " + emp::to_string(i), {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
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
      }, 1, "Terminal value " + emp::to_string(i), {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  }
}

// Problem setups
void ProgramSynthesisExperiment::SetupProblem_NumberIO() {
  std::cout << "Setting up problem: NumberIO." << std::endl;

  using prob_input_t = typename ProblemUtilities_NumberIO::input_t;
  using prob_output_t = typename ProblemUtilities_NumberIO::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  USES_VECTOR_INSTRUCTIONS = false;
  USES_STRING_INSTRUCTIONS = false;

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
  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      
      inst_lib->AddInst("LoadInt-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadInt_NumberIO__TAG_ARGS(hw, inst); }, 1, "LoadInt-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadDouble-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadDouble_NumberIO__TAG_ARGS(hw, inst); }, 1, "LoadDouble-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_NumberIO__TAG_ARGS(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_NumberIO__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      
      inst_lib->AddInst("LoadInt-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadInt_NumberIO__NUM_ARGS(hw, inst); }, 1, "LoadInt-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadDouble-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadDouble_NumberIO__NUM_ARGS(hw, inst); }, 1, "LoadDouble-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_NumberIO__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_NumberIO__NUM_ARGS(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      // Tag-based arguments
      inst_lib->AddInst("LoadInt-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadInt_NumberIO__TAG_ARGS(hw, inst); }, 1, "LoadInt-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadDouble-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadDouble_NumberIO__TAG_ARGS(hw, inst); }, 1, "LoadDouble-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_NumberIO__TAG_ARGS(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_NumberIO__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      // Numeric arguments
      inst_lib->AddInst("LoadInt-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadInt_NumberIO__NUM_ARGS(hw, inst); }, 1, "LoadInt-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadDouble-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadDouble_NumberIO__NUM_ARGS(hw, inst); }, 1, "LoadDouble-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_NumberIO__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_NumberIO__NUM_ARGS(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      break;
    }
  }
}

void ProgramSynthesisExperiment::SetupProblem_SmallOrLarge() {
  std::cout << "Setting up problem: SmallOrLarge." << std::endl;

  using prob_input_t = typename ProblemUtilities_SmallOrLarge::input_t;
  using prob_output_t = typename ProblemUtilities_SmallOrLarge::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  USES_VECTOR_INSTRUCTIONS = false;
  USES_STRING_INSTRUCTIONS = false;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_SmallOrLarge.GetTrainingSet().LoadTestCases(training_examples_fpath);
  prob_utils_SmallOrLarge.GetTestingSet().LoadTestCases(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_SmallOrLarge.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_SmallOrLarge.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_SmallOrLarge.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_SmallOrLarge.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_SmallOrLarge.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input);
    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_SmallOrLarge.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_SmallOrLarge.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_SmallOrLarge.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_SmallOrLarge.CalcScorePassFail(correct_output, prob_utils_SmallOrLarge.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);
  
  inst_lib->AddInst("SubmitSmall", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitSmall_SmallOrLarge(hw, inst); }, 0, "SubmitSmall", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SubmitLarge", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitLarge_SmallOrLarge(hw, inst); }, 0, "SubmitLarge", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SubmitNone", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNone_SmallOrLarge(hw, inst); }, 0, "SubmitNone", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  
  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      inst_lib->AddInst("LoadInt-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadInt_SmallOrLarge__TAG_ARGS(hw, inst); }, 1, "LoadInt-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      inst_lib->AddInst("LoadInt-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadInt_SmallOrLarge__NUM_ARGS(hw, inst); }, 1, "LoadInt-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      inst_lib->AddInst("LoadInt-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadInt_SmallOrLarge__TAG_ARGS(hw, inst); }, 1, "LoadInt-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadInt-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadInt_SmallOrLarge__NUM_ARGS(hw, inst); }, 1, "LoadInt-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      break;
    }
  }
}

void ProgramSynthesisExperiment::SetupProblem_ForLoopIndex() {
  std::cout << "Setting up problem: ForLoopIndex." << std::endl;

  using prob_input_t = typename ProblemUtilities_ForLoopIndex::input_t;
  using prob_output_t = typename ProblemUtilities_ForLoopIndex::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  USES_VECTOR_INSTRUCTIONS = false;
  USES_STRING_INSTRUCTIONS = false;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_ForLoopIndex.GetTrainingSet().LoadTestCases(training_examples_fpath);
  prob_utils_ForLoopIndex.GetTestingSet().LoadTestCases(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_ForLoopIndex.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_ForLoopIndex.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_ForLoopIndex.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 3);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_ForLoopIndex.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_ForLoopIndex.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
      wmem.Set(2, input[2]);
    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_ForLoopIndex.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_ForLoopIndex.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_ForLoopIndex.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_ForLoopIndex.CalcScoreGradient(correct_output, prob_utils_ForLoopIndex.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);

  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      inst_lib->AddInst("LoadStart-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStart_ForLoopIndex__TAG_ARGS(hw, inst); }, 1, "LoadStart-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadEnd-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadEnd_ForLoopIndex__TAG_ARGS(hw, inst); }, 1, "LoadEnd-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStep-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStep_ForLoopIndex__TAG_ARGS(hw, inst); }, 1, "LoadStep-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_ForLoopIndex__TAG_ARGS(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_ForLoopIndex__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
       
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      inst_lib->AddInst("LoadStart-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStart_ForLoopIndex__NUM_ARGS(hw, inst); }, 1, "LoadStart-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadEnd-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadEnd_ForLoopIndex__NUM_ARGS(hw, inst); }, 1, "LoadEnd-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStep-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStep_ForLoopIndex__NUM_ARGS(hw, inst); }, 1, "LoadStep-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_ForLoopIndex__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_ForLoopIndex__NUM_ARGS(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      // tag-based arguments
      inst_lib->AddInst("LoadStart-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStart_ForLoopIndex__TAG_ARGS(hw, inst); }, 1, "LoadStart-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadEnd-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadEnd_ForLoopIndex__TAG_ARGS(hw, inst); }, 1, "LoadEnd-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStep-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStep_ForLoopIndex__TAG_ARGS(hw, inst); }, 1, "LoadStep-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_ForLoopIndex__TAG_ARGS(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_ForLoopIndex__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      // Numeric arguments
      inst_lib->AddInst("LoadStart-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStart_ForLoopIndex__NUM_ARGS(hw, inst); }, 1, "LoadStart-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadEnd-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadEnd_ForLoopIndex__NUM_ARGS(hw, inst); }, 1, "LoadEnd-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStep-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStep_ForLoopIndex__NUM_ARGS(hw, inst); }, 1, "LoadStep-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_ForLoopIndex__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_ForLoopIndex__NUM_ARGS(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }   

      break;
    }
  }
}

void ProgramSynthesisExperiment::SetupProblem_CompareStringLengths() {
  std::cout << "Setting up problem: CompareStringLengths." << std::endl;

  using prob_input_t = typename ProblemUtilities_CompareStringLengths::input_t;
  using prob_output_t = typename ProblemUtilities_CompareStringLengths::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  USES_VECTOR_INSTRUCTIONS = false;
  USES_STRING_INSTRUCTIONS = true;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_CompareStringLengths.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  prob_utils_CompareStringLengths.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_CompareStringLengths.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_CompareStringLengths.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_CompareStringLengths.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 3);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_CompareStringLengths.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_CompareStringLengths.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
      wmem.Set(2, input[2]);
    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_CompareStringLengths.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_CompareStringLengths.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_CompareStringLengths.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_CompareStringLengths.CalcScorePassFail(correct_output, prob_utils_CompareStringLengths.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);

  inst_lib->AddInst("SubmitTrue", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitTrue_CompareStringLengths(hw, inst); }, 0, "SubmitTrue", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SubmitFalse", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitFalse_CompareStringLengths(hw, inst); }, 0, "SubmitFalse", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {

      inst_lib->AddInst("LoadStr1-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr1_CompareStringLengths__TAG_ARGS(hw, inst); }, 1, "LoadStr1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStr2-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr2_CompareStringLengths__TAG_ARGS(hw, inst); }, 1, "LoadStr2", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStr3-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr3_CompareStringLengths__TAG_ARGS(hw, inst); }, 1, "LoadStr3", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_CompareStringLengths__TAG_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_CompareStringLengths__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {

      inst_lib->AddInst("LoadStr1-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr1_CompareStringLengths__NUM_ARGS(hw, inst); }, 1, "LoadStr1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStr2-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr2_CompareStringLengths__NUM_ARGS(hw, inst); }, 1, "LoadStr2", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStr3-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr3_CompareStringLengths__NUM_ARGS(hw, inst); }, 1, "LoadStr3", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_CompareStringLengths__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_CompareStringLengths__NUM_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      // Tag-based arguments
      inst_lib->AddInst("LoadStr1-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr1_CompareStringLengths__TAG_ARGS(hw, inst); }, 1, "LoadStr1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStr2-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr2_CompareStringLengths__TAG_ARGS(hw, inst); }, 1, "LoadStr2", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStr3-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr3_CompareStringLengths__TAG_ARGS(hw, inst); }, 1, "LoadStr3", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_CompareStringLengths__TAG_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_CompareStringLengths__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      // Numeric arguments
      inst_lib->AddInst("LoadStr1-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr1_CompareStringLengths__NUM_ARGS(hw, inst); }, 1, "LoadStr1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStr2-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr2_CompareStringLengths__NUM_ARGS(hw, inst); }, 1, "LoadStr2", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadStr3-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStr3_CompareStringLengths__NUM_ARGS(hw, inst); }, 1, "LoadStr3", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_CompareStringLengths__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_CompareStringLengths__NUM_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      break;
    }
  }
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
  std::cout << "Setting up problem: StringLengthsBackwards." << std::endl;

  using prob_input_t = typename ProblemUtilities_StringLengthsBackwards::input_t;
  using prob_output_t = typename ProblemUtilities_StringLengthsBackwards::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  USES_VECTOR_INSTRUCTIONS = true;
  USES_STRING_INSTRUCTIONS = true;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_StringLengthsBackwards.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  prob_utils_StringLengthsBackwards.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_StringLengthsBackwards.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_StringLengthsBackwards.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_StringLengthsBackwards.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 3);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_StringLengthsBackwards.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_StringLengthsBackwards.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
      wmem.Set(2, input[2]);
    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_StringLengthsBackwards.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_StringLengthsBackwards.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_StringLengthsBackwards.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_StringLengthsBackwards.CalcScorePassFail(correct_output, prob_utils_StringLengthsBackwards.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);

  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      
      inst_lib->AddInst("LoadStrVec-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStrVec_StringLengthsBackwards__TAG_ARGS(hw, inst); }, 1, "LoadStrVec", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_StringLengthsBackwards__TAG_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_StringLengthsBackwards__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      
      inst_lib->AddInst("LoadStrVec-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStrVec_StringLengthsBackwards__NUM_ARGS(hw, inst); }, 1, "LoadStrVec", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_StringLengthsBackwards__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_StringLengthsBackwards__NUM_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      // Tag-based arguments
      inst_lib->AddInst("LoadStrVec-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStrVec_StringLengthsBackwards__TAG_ARGS(hw, inst); }, 1, "LoadStrVec", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_StringLengthsBackwards__TAG_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_StringLengthsBackwards__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      // Numeric arguments
      inst_lib->AddInst("LoadStrVec-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadStrVec_StringLengthsBackwards__NUM_ARGS(hw, inst); }, 1, "LoadStrVec", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_StringLengthsBackwards__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_StringLengthsBackwards__NUM_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      break;
    }
  }
}

void ProgramSynthesisExperiment::SetupProblem_LastIndexOfZero() {
  std::cout << "Setting up problem: LastIndexOfZero." << std::endl;

  using prob_input_t = typename ProblemUtilities_LastIndexOfZero::input_t;
  using prob_output_t = typename ProblemUtilities_LastIndexOfZero::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  USES_VECTOR_INSTRUCTIONS = true;
  USES_STRING_INSTRUCTIONS = false;

  prob_utils_LastIndexOfZero.MAX_ERROR = PROB_LAST_INDEX_OF_ZERO__MAX_VEC_LEN;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_LastIndexOfZero.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  prob_utils_LastIndexOfZero.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_LastIndexOfZero.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_LastIndexOfZero.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_LastIndexOfZero.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_LastIndexOfZero.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_LastIndexOfZero.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input);
    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_LastIndexOfZero.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_LastIndexOfZero.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_LastIndexOfZero.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_LastIndexOfZero.CalcScoreGradient(correct_output, prob_utils_LastIndexOfZero.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);

  std::cout << "Problem-specific instructions not yet implemented. Exiting." << std::endl;
  exit(-1);
  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      inst_lib->AddInst("LoadVec-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec_LastIndexOfZero__TAG_ARGS(hw, inst); }, 1, "LoadVec", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_LastIndexOfZero__TAG_ARGS(hw, inst); }, 1, "SubmitNum", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_LastIndexOfZero__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      inst_lib->AddInst("LoadVec-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec_LastIndexOfZero__NUM_ARGS(hw, inst); }, 1, "LoadVec", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_LastIndexOfZero__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_LastIndexOfZero__NUM_ARGS(hw, inst); }, 1, "SubmitNum", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      // Tag-based arguments
      inst_lib->AddInst("LoadVec-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec_LastIndexOfZero__TAG_ARGS(hw, inst); }, 1, "LoadVec", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_LastIndexOfZero__TAG_ARGS(hw, inst); }, 1, "SubmitNum", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_LastIndexOfZero__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      // Numeric arguments
      inst_lib->AddInst("LoadVec-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec_LastIndexOfZero__NUM_ARGS(hw, inst); }, 1, "LoadVec", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_LastIndexOfZero__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_LastIndexOfZero__NUM_ARGS(hw, inst); }, 1, "SubmitNum", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      break;
    }
  }
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
  std::cout << "Setting up problem: MirrorImage." << std::endl;

  using prob_input_t = typename ProblemUtilities_MirrorImage::input_t;
  using prob_output_t = typename ProblemUtilities_MirrorImage::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_MirrorImage.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  prob_utils_MirrorImage.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_MirrorImage.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_MirrorImage.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_MirrorImage.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_MirrorImage.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_MirrorImage.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_MirrorImage.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_MirrorImage.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_MirrorImage.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_MirrorImage.CalcScorePassFail(correct_output, prob_utils_MirrorImage.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);

  inst_lib->AddInst("SubmitTrue", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitTrue_MirrorImage(hw, inst); }, 0, "SubmitTrue", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SubmitFalse", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitFalse_MirrorImage(hw, inst); }, 0, "SubmitFalse", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  
  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      inst_lib->AddInst("LoadVec1-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec1_MirrorImage__TAG_ARGS(hw, inst); }, 1, "LoadVec1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadVec2-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec2_MirrorImage__TAG_ARGS(hw, inst); }, 1, "LoadVec2", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_MirrorImage__TAG_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_MirrorImage__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      inst_lib->AddInst("LoadVec1-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec1_MirrorImage__NUM_ARGS(hw, inst); }, 1, "LoadVec1", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadVec2-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec2_MirrorImage__NUM_ARGS(hw, inst); }, 1, "LoadVec2", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_MirrorImage__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_MirrorImage__NUM_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      // Tag-based arguments
      inst_lib->AddInst("LoadVec1-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec1_MirrorImage__TAG_ARGS(hw, inst); }, 1, "LoadVec1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadVec2-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec2_MirrorImage__TAG_ARGS(hw, inst); }, 1, "LoadVec2", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_MirrorImage__TAG_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_MirrorImage__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      // Numeric arguments
      inst_lib->AddInst("LoadVec1-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec1_MirrorImage__NUM_ARGS(hw, inst); }, 1, "LoadVec1", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadVec2-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec2_MirrorImage__NUM_ARGS(hw, inst); }, 1, "LoadVec2", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_MirrorImage__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitVal", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVal-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVal_MirrorImage__NUM_ARGS(hw, inst); }, 1, "SubmitVal", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }      
      break;
    }
  }
}

void ProgramSynthesisExperiment::SetupProblem_SuperAnagrams() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

void ProgramSynthesisExperiment::SetupProblem_SumOfSquares() {
  std::cout << "Setting up problem: SumOfSquares." << std::endl;

  using prob_input_t = typename ProblemUtilities_SumOfSquares::input_t;
  using prob_output_t = typename ProblemUtilities_SumOfSquares::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  USES_VECTOR_INSTRUCTIONS = false;
  USES_STRING_INSTRUCTIONS = false;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_SumOfSquares.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  prob_utils_SumOfSquares.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_SumOfSquares.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_SumOfSquares.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_SumOfSquares.ResetTestEval();
    prob_utils_SumOfSquares.MAX_ERROR = (int)(prob_utils_SumOfSquares.GetTestOutput(eval_util.use_training_set, eval_util.current_testID) * 0.5);
    
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_SumOfSquares.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_SumOfSquares.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input);
    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_SumOfSquares.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_SumOfSquares.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_SumOfSquares.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_SumOfSquares.CalcScoreGradient(correct_output, prob_utils_SumOfSquares.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);

  // todo
  std::cout << "Problem-specific instructions not yet implemented. Exiting." << std::endl;
  exit(-1);
  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      break;
    }
  }
}

void ProgramSynthesisExperiment::SetupProblem_VectorsSummed() {
  std::cout << "Setting up problem: VectorsSummed." << std::endl;

  using prob_input_t = typename ProblemUtilities_VectorsSummed::input_t;
  using prob_output_t = typename ProblemUtilities_VectorsSummed::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  USES_VECTOR_INSTRUCTIONS = true;
  USES_STRING_INSTRUCTIONS = false;

  prob_utils_VectorsSummed.MAX_NUM = PROB_VECTORS_SUMMED__MAX_NUM;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_VectorsSummed.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  prob_utils_VectorsSummed.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_VectorsSummed.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_VectorsSummed.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_VectorsSummed.ResetTestEval();
    prob_utils_VectorsSummed.MAX_ERROR = (2*PROB_VECTORS_SUMMED__MAX_NUM) * prob_utils_VectorsSummed.GetTestOutput(eval_util.use_training_set, eval_util.current_testID).size();
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_VectorsSummed.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_VectorsSummed.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_VectorsSummed.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_VectorsSummed.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_VectorsSummed.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_VectorsSummed.CalcScoreGradient(correct_output, prob_utils_VectorsSummed.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);

  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      inst_lib->AddInst("LoadVec1-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec1_VectorsSummed__TAG_ARGS(hw, inst); }, 1, "LoadVec1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadVec2-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec2_VectorsSummed__TAG_ARGS(hw, inst); }, 1, "LoadVec2", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVec-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVec_VectorsSummed__TAG_ARGS(hw, inst); }, 1, "SubmitVec", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVec-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVec_VectorsSummed__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitVec", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      inst_lib->AddInst("LoadVec1-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec1_VectorsSummed__NUM_ARGS(hw, inst); }, 1, "LoadVec1", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadVec2-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec2_VectorsSummed__NUM_ARGS(hw, inst); }, 1, "LoadVec2", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVec-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVec_VectorsSummed__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitVec", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING}); 
      } else {
        inst_lib->AddInst("SubmitVec-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVec_VectorsSummed__NUM_ARGS(hw, inst); }, 1, "SubmitVec", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      // Tag-based arguments
      inst_lib->AddInst("LoadVec1-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec1_VectorsSummed__TAG_ARGS(hw, inst); }, 1, "LoadVec1", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadVec2-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec2_VectorsSummed__TAG_ARGS(hw, inst); }, 1, "LoadVec2", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVec-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVec_VectorsSummed__TAG_ARGS(hw, inst); }, 1, "SubmitVec", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitVec-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVec_VectorsSummed__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitVec", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      // Numeric arguments
      inst_lib->AddInst("LoadVec1-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec1_VectorsSummed__NUM_ARGS(hw, inst); }, 1, "LoadVec1", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadVec2-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadVec2_VectorsSummed__NUM_ARGS(hw, inst); }, 1, "LoadVec2", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitVec-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVec_VectorsSummed__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitVec", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING}); 
      } else {
        inst_lib->AddInst("SubmitVec-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitVec_VectorsSummed__NUM_ARGS(hw, inst); }, 1, "SubmitVec", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }      
      break;
    }
  }
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
  std::cout << "Setting up problem: Grade." << std::endl;

  using prob_input_t = typename ProblemUtilities_Grade::input_t;
  using prob_output_t = typename ProblemUtilities_Grade::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  USES_VECTOR_INSTRUCTIONS = false;
  USES_STRING_INSTRUCTIONS = false;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_Grade.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  prob_utils_Grade.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_Grade.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_Grade.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_Grade.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_Grade.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_Grade.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
      wmem.Set(2, input[2]);
      wmem.Set(3, input[3]);
      wmem.Set(4, input[4]);
    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_Grade.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_Grade.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_Grade.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_Grade.CalcScorePassFail(correct_output, prob_utils_Grade.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);

  inst_lib->AddInst("SubmitA", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitA_Grade(hw, inst); }, 0, "SubmitA", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SubmitB", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitB_Grade(hw, inst); }, 0, "SubmitB", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SubmitC", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitC_Grade(hw, inst); }, 0, "SubmitC", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SubmitD", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitD_Grade(hw, inst); }, 0, "SubmitD", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
  inst_lib->AddInst("SubmitF", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitF_Grade(hw, inst); }, 0, "SubmitF", {inst_prop_t::NO_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      inst_lib->AddInst("LoadThreshA-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshA_Grade__TAG_ARGS(hw, inst); }, 1, "LoadThreshA-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshB-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshB_Grade__TAG_ARGS(hw, inst); }, 1, "LoadThreshB-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshC-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshC_Grade__TAG_ARGS(hw, inst); }, 1, "LoadThreshC-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshD-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshD_Grade__TAG_ARGS(hw, inst); }, 1, "LoadThreshD-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadGrade-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadGrade_Grade__TAG_ARGS(hw, inst); }, 1, "LoadGrade-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      inst_lib->AddInst("LoadThreshA-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshA_Grade__NUM_ARGS(hw, inst); }, 1, "LoadThreshA-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshB-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshB_Grade__NUM_ARGS(hw, inst); }, 1, "LoadThreshB-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshC-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshC_Grade__NUM_ARGS(hw, inst); }, 1, "LoadThreshC-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshD-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshD_Grade__NUM_ARGS(hw, inst); }, 1, "LoadThreshD-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadGrade-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadGrade_Grade__NUM_ARGS(hw, inst); }, 1, "LoadGrade-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      inst_lib->AddInst("LoadThreshA-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshA_Grade__TAG_ARGS(hw, inst); }, 1, "LoadThreshA-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshB-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshB_Grade__TAG_ARGS(hw, inst); }, 1, "LoadThreshB-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshC-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshC_Grade__TAG_ARGS(hw, inst); }, 1, "LoadThreshC-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshD-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshD_Grade__TAG_ARGS(hw, inst); }, 1, "LoadThreshD-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadGrade-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadGrade_Grade__TAG_ARGS(hw, inst); }, 1, "LoadGrade-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 

      inst_lib->AddInst("LoadThreshA-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshA_Grade__NUM_ARGS(hw, inst); }, 1, "LoadThreshA-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshB-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshB_Grade__NUM_ARGS(hw, inst); }, 1, "LoadThreshB-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshC-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshC_Grade__NUM_ARGS(hw, inst); }, 1, "LoadThreshC-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadThreshD-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadThreshD_Grade__NUM_ARGS(hw, inst); }, 1, "LoadThreshD-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      inst_lib->AddInst("LoadGrade-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadGrade_Grade__NUM_ARGS(hw, inst); }, 1, "LoadGrade-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC}); 
      break;
    }
  }
}

void ProgramSynthesisExperiment::SetupProblem_Median() {
  std::cout << "Setting up problem: Median." << std::endl;

  using prob_input_t = typename ProblemUtilities_Median::input_t;
  using prob_output_t = typename ProblemUtilities_Median::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  USES_VECTOR_INSTRUCTIONS = false;
  USES_STRING_INSTRUCTIONS = false;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_Median.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  prob_utils_Median.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_Median.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_Median.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_Median.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_Median.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_Median.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
      wmem.Set(2, input[2]);
    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_Median.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_Median.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_Median.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_Median.CalcScorePassFail(correct_output, prob_utils_Median.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);
  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      inst_lib->AddInst("LoadNum1-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum1_Median__TAG_ARGS(hw, inst); }, 1, "LoadNum1-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum2-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum2_Median__TAG_ARGS(hw, inst); }, 1, "LoadNum2-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum3-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum3_Median__TAG_ARGS(hw, inst); }, 1, "LoadNum3-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Median__TAG_ARGS(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Median__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      inst_lib->AddInst("LoadNum1-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum1_Median__NUM_ARGS(hw, inst); }, 1, "LoadNum1-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum2-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum2_Median__NUM_ARGS(hw, inst); }, 1, "LoadNum2-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum3-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum3_Median__NUM_ARGS(hw, inst); }, 1, "LoadNum3-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Median__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Median__NUM_ARGS(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      
      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      // Tag-based arguments
      inst_lib->AddInst("LoadNum1-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum1_Median__TAG_ARGS(hw, inst); }, 1, "LoadNum1-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum2-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum2_Median__TAG_ARGS(hw, inst); }, 1, "LoadNum2-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum3-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum3_Median__TAG_ARGS(hw, inst); }, 1, "LoadNum3-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Median__TAG_ARGS(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Median__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      
      // Numeric arguments
      inst_lib->AddInst("LoadNum1-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum1_Median__NUM_ARGS(hw, inst); }, 1, "LoadNum1-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum2-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum2_Median__NUM_ARGS(hw, inst); }, 1, "LoadNum2-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum3-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum3_Median__NUM_ARGS(hw, inst); }, 1, "LoadNum3-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Median__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Median__NUM_ARGS(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

     break;
    }
  }
}

void ProgramSynthesisExperiment::SetupProblem_Smallest() {
  std::cout << "Setting up problem: Smallest." << std::endl;

  using prob_input_t = typename ProblemUtilities_Smallest::input_t;
  using prob_output_t = typename ProblemUtilities_Smallest::output_t;
  using testcase_set_t = TestCaseSet<prob_input_t,prob_output_t>;

  USES_VECTOR_INSTRUCTIONS = false;
  USES_STRING_INSTRUCTIONS = false;
 
  // Load training and testing examples from file.
  if (BENCHMARK_DATA_DIR.back() != '/') BENCHMARK_DATA_DIR += '/';
  std::string training_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTrainingSetFilename();  
  std::string testing_examples_fpath = BENCHMARK_DATA_DIR + problems.at(PROBLEM).GetTestingSetFilename();  
  prob_utils_Smallest.GetTrainingSet().LoadTestCasesWithCSVReader(training_examples_fpath);
  prob_utils_Smallest.GetTestingSet().LoadTestCasesWithCSVReader(testing_examples_fpath);
  TRAINING_SET_SIZE = prob_utils_Smallest.GetTrainingSet().GetSize();
  TESTING_SET_SIZE = prob_utils_Smallest.GetTestingSet().GetSize();
  std::cout << "Loaded TRAINING set size = " << TRAINING_SET_SIZE << std::endl;
  std::cout << "Loaded TESTING set size  = " << TESTING_SET_SIZE << std::endl;

  // Tell experiment how to configure hardware inputs when running a program against a test.
  begin_program_test.AddAction([this](prog_org_t & prog_org) {
    // Reset evaluation utilities.
    prob_utils_Smallest.ResetTestEval();
    emp_assert(eval_hardware->GetMemSize() >= 1);
    // Configure inputs.
    if (eval_hardware->GetCallStackSize()) {
      
      // Are we using the training set or testing set?
      emp::Ptr<testcase_set_t> test_set_ptr;
      if (eval_util.use_training_set) test_set_ptr = &prob_utils_Smallest.training_set; // todo - confirm this is okay
      else test_set_ptr = &prob_utils_Smallest.testing_set;

      emp_assert(eval_util.current_testID < test_set_ptr->GetSize());
      prob_input_t & input = test_set_ptr->GetInput(eval_util.current_testID);
      hardware_t::CallState & state = eval_hardware->GetCurCallState();
      hardware_t::Memory & wmem = state.GetWorkingMem();

      // Set hardware inputs.
      wmem.Set(0, input[0]);
      wmem.Set(1, input[1]);
      wmem.Set(2, input[2]);
      wmem.Set(3, input[3]);
    }
  });

  // Tell the experiment how to calculate test results.
  CalcProgramResultOnTest = [this](prog_org_t & prog_org) {
    // Are we using the training set or testing set?
    emp::Ptr<testcase_set_t> test_set_ptr;
    if (eval_util.use_training_set) test_set_ptr = &prob_utils_Smallest.training_set; // todo - confirm this is okay
    else test_set_ptr = &prob_utils_Smallest.testing_set;

    prob_output_t & correct_output = test_set_ptr->GetOutput(eval_util.current_testID);
    
    TestResult result;
    if (!prob_utils_Smallest.submitted) {
      result.score = 0; result.pass = false; result.sub = false;
    } else {
      std::pair<double, bool> r(prob_utils_Smallest.CalcScorePassFail(correct_output, prob_utils_Smallest.submitted_val));
      result.score = r.first; result.pass = r.second; result.sub = true;
    }
    return result;
  };

  // Add problem-specific instructions. (Terminals)
  AddNumericTerminals(0, 16);

  switch (PROGRAM_ARGUMENT_MODE) {
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::TAG_ONLY: {
      inst_lib->AddInst("LoadNum1-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum1_Smallest__TAG_ARGS(hw, inst); }, 1, "LoadNum1-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum2-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum2_Smallest__TAG_ARGS(hw, inst); }, 1, "LoadNum2-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum3-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum3_Smallest__TAG_ARGS(hw, inst); }, 1, "LoadNum3-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum4-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum4_Smallest__TAG_ARGS(hw, inst); }, 1, "LoadNum4-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Smallest__TAG_ARGS(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Smallest__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::NUMERIC_ONLY: {
      inst_lib->AddInst("LoadNum1-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum1_Smallest__NUM_ARGS(hw, inst); }, 1, "LoadNum1-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum2-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum2_Smallest__NUM_ARGS(hw, inst); }, 1, "LoadNum2-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum3-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum3_Smallest__NUM_ARGS(hw, inst); }, 1, "LoadNum3-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum4-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum4_Smallest__NUM_ARGS(hw, inst); }, 1, "LoadNum4-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Smallest__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Smallest__NUM_ARGS(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      break;
    }
    case (size_t)PROGRAM_ARGUMENT_MODE_TYPE::BOTH: {
      // tag-based arguments
      inst_lib->AddInst("LoadNum1-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum1_Smallest__TAG_ARGS(hw, inst); }, 1, "LoadNum1-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum2-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum2_Smallest__TAG_ARGS(hw, inst); }, 1, "LoadNum2-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum3-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum3_Smallest__TAG_ARGS(hw, inst); }, 1, "LoadNum3-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum4-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum4_Smallest__TAG_ARGS(hw, inst); }, 1, "LoadNum4-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      
      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Smallest__TAG_ARGS(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Tag", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Smallest__TAG_ARGS_NO_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Tag", {inst_prop_t::TAG_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }

      // Numeric arguments
      inst_lib->AddInst("LoadNum1-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum1_Smallest__NUM_ARGS(hw, inst); }, 1, "LoadNum1-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum2-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum2_Smallest__NUM_ARGS(hw, inst); }, 1, "LoadNum2-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum3-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum3_Smallest__NUM_ARGS(hw, inst); }, 1, "LoadNum3-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});
      inst_lib->AddInst("LoadNum4-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_LoadNum4_Smallest__NUM_ARGS(hw, inst); }, 1, "LoadNum4-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_AGNOSTIC});

      if (PROGRAM_ARGUMENTS_TYPE_SEARCH) {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Smallest__NUM_ARGS_WITH_TYPE_SEARCH(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_SEARCHING});
      } else {
        inst_lib->AddInst("SubmitNum-Num", [this](hardware_t & hw, const inst_t & inst) { this->Inst_SubmitNum_Smallest__NUM_ARGS(hw, inst); }, 1, "SubmitNum-Num", {inst_prop_t::NUM_ARGS, inst_prop_t::MEM_TYPE_NO_SEARCHING});
      }
      break;
    }
  }
}

void ProgramSynthesisExperiment::SetupProblem_Syllables() {
  std::cout << "Problem setup is not yet implemented. Exiting." << std::endl;
  exit(-1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Instruction implementations /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Number-io
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ProgramSynthesisExperiment::Inst_LoadInt_NumberIO__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_NumberIO.GetTestInput(eval_util.use_training_set, eval_util.current_testID).first);
}

void ProgramSynthesisExperiment::Inst_LoadDouble_NumberIO__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_NumberIO.GetTestInput(eval_util.use_training_set, eval_util.current_testID).second);
}

void ProgramSynthesisExperiment::Inst_SubmitNum_NumberIO__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_NumberIO.Submit(wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_LoadInt_NumberIO__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_NumberIO.GetTestInput(eval_util.use_training_set, eval_util.current_testID).first);
}

void ProgramSynthesisExperiment::Inst_LoadDouble_NumberIO__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_NumberIO.GetTestInput(eval_util.use_training_set, eval_util.current_testID).second);
}

void ProgramSynthesisExperiment::Inst_SubmitNum_NumberIO__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_NumberIO.Submit(wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_NumberIO__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_NumberIO.Submit(wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_NumberIO__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_nums[0], hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_NumberIO.Submit(wmem.AccessVal(posA).GetNum());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Small or large
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ProgramSynthesisExperiment::Inst_LoadInt_SmallOrLarge__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_SmallOrLarge.GetTestInput(eval_util.use_training_set, eval_util.current_testID));
}

void ProgramSynthesisExperiment::Inst_LoadInt_SmallOrLarge__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_SmallOrLarge.GetTestInput(eval_util.use_training_set, eval_util.current_testID));
}

void ProgramSynthesisExperiment::Inst_SubmitSmall_SmallOrLarge(hardware_t & hw, const inst_t & inst) { 
  prob_utils_SmallOrLarge.Submit("small");
}

void ProgramSynthesisExperiment::Inst_SubmitLarge_SmallOrLarge(hardware_t & hw, const inst_t & inst) { 
  prob_utils_SmallOrLarge.Submit("large");
}

void ProgramSynthesisExperiment::Inst_SubmitNone_SmallOrLarge(hardware_t & hw, const inst_t & inst) { 
  prob_utils_SmallOrLarge.Submit("");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// For Loop Index
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// -- Tag-based arguments --
void ProgramSynthesisExperiment::Inst_LoadStart_ForLoopIndex__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_ForLoopIndex.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadEnd_ForLoopIndex__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_ForLoopIndex.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}

void ProgramSynthesisExperiment::Inst_LoadStep_ForLoopIndex__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_ForLoopIndex.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[2]);
}

void ProgramSynthesisExperiment::Inst_SubmitNum_ForLoopIndex__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_ForLoopIndex.Submit((int)wmem.AccessVal(posA).GetNum());
}

// -- Numeric arguments --
void ProgramSynthesisExperiment::Inst_LoadStart_ForLoopIndex__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_ForLoopIndex.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadEnd_ForLoopIndex__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_ForLoopIndex.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]); 
}

void ProgramSynthesisExperiment::Inst_LoadStep_ForLoopIndex__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_ForLoopIndex.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[2]); 
}

void ProgramSynthesisExperiment::Inst_SubmitNum_ForLoopIndex__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_ForLoopIndex.Submit((int)wmem.AccessVal(posA).GetNum()); 
}

void ProgramSynthesisExperiment::Inst_SubmitNum_ForLoopIndex__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_ForLoopIndex.Submit((int)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_ForLoopIndex__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_nums[0], hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_ForLoopIndex.Submit((int)wmem.AccessVal(posA).GetNum());
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Grade
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// -- Tag-based arguments --
void ProgramSynthesisExperiment::Inst_LoadThreshA_Grade__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Grade.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadThreshB_Grade__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Grade.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}

void ProgramSynthesisExperiment::Inst_LoadThreshC_Grade__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Grade.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[2]);
}

void ProgramSynthesisExperiment::Inst_LoadThreshD_Grade__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Grade.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[3]);
}

void ProgramSynthesisExperiment::Inst_LoadGrade_Grade__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Grade.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[4]);
}

// -- Numeric arguments --
void ProgramSynthesisExperiment::Inst_LoadThreshA_Grade__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Grade.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadThreshB_Grade__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Grade.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}

void ProgramSynthesisExperiment::Inst_LoadThreshC_Grade__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Grade.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[2]);
}

void ProgramSynthesisExperiment::Inst_LoadThreshD_Grade__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Grade.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[3]);
}

void ProgramSynthesisExperiment::Inst_LoadGrade_Grade__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Grade.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[4]);
}

void ProgramSynthesisExperiment::Inst_SubmitA_Grade(hardware_t & hw, const inst_t & inst) { prob_utils_Grade.Submit("A"); }
void ProgramSynthesisExperiment::Inst_SubmitB_Grade(hardware_t & hw, const inst_t & inst) { prob_utils_Grade.Submit("B"); }
void ProgramSynthesisExperiment::Inst_SubmitC_Grade(hardware_t & hw, const inst_t & inst) { prob_utils_Grade.Submit("C"); }
void ProgramSynthesisExperiment::Inst_SubmitD_Grade(hardware_t & hw, const inst_t & inst) { prob_utils_Grade.Submit("D"); }
void ProgramSynthesisExperiment::Inst_SubmitF_Grade(hardware_t & hw, const inst_t & inst) { prob_utils_Grade.Submit("F"); }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Median
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tag-based arguments
void ProgramSynthesisExperiment::Inst_LoadNum1_Median__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Median.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadNum2_Median__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Median.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}

void ProgramSynthesisExperiment::Inst_LoadNum3_Median__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Median.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[2]);
}

void ProgramSynthesisExperiment::Inst_SubmitNum_Median__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_Median.Submit((int)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_Median__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_Median.Submit((int)wmem.AccessVal(posA).GetNum());
}

// -- Numeric arguments --
void ProgramSynthesisExperiment::Inst_LoadNum1_Median__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Median.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadNum2_Median__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Median.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}

void ProgramSynthesisExperiment::Inst_LoadNum3_Median__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Median.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[2]);
}

void ProgramSynthesisExperiment::Inst_SubmitNum_Median__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_Median.Submit((int)wmem.AccessVal(posA).GetNum());   
}

void ProgramSynthesisExperiment::Inst_SubmitNum_Median__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_nums[0], hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_Median.Submit((int)wmem.AccessVal(posA).GetNum());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Smallest
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ProgramSynthesisExperiment::Inst_LoadNum1_Smallest__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Smallest.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);  
}
void ProgramSynthesisExperiment::Inst_LoadNum2_Smallest__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Smallest.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);  
}
void ProgramSynthesisExperiment::Inst_LoadNum3_Smallest__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Smallest.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[2]);  
}
void ProgramSynthesisExperiment::Inst_LoadNum4_Smallest__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Smallest.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[3]);  
}
void ProgramSynthesisExperiment::Inst_SubmitNum_Smallest__TAG_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_Smallest.Submit((int)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_Smallest__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_Smallest.Submit((int)wmem.AccessVal(posA).GetNum()); 
}

// -- Numeric arguments --
void ProgramSynthesisExperiment::Inst_LoadNum1_Smallest__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Smallest.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}
void ProgramSynthesisExperiment::Inst_LoadNum2_Smallest__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Smallest.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}
void ProgramSynthesisExperiment::Inst_LoadNum3_Smallest__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Smallest.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[2]);
}
void ProgramSynthesisExperiment::Inst_LoadNum4_Smallest__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_Smallest.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[3]);  
}
void ProgramSynthesisExperiment::Inst_SubmitNum_Smallest__NUM_ARGS(hardware_t & hw, const inst_t & inst) { 
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_Smallest.Submit((int)wmem.AccessVal(posA).GetNum());  
}

void ProgramSynthesisExperiment::Inst_SubmitNum_Smallest__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_nums[0], hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_Smallest.Submit((int)wmem.AccessVal(posA).GetNum());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// -- Compare string lengths --
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ProgramSynthesisExperiment::Inst_SubmitTrue_CompareStringLengths(hardware_t & hw, const inst_t & inst) {
  prob_utils_CompareStringLengths.Submit(true);
}

void ProgramSynthesisExperiment::Inst_SubmitFalse_CompareStringLengths(hardware_t & hw, const inst_t & inst) {
  prob_utils_CompareStringLengths.Submit(false);
}

// Tag-based args
void ProgramSynthesisExperiment::Inst_LoadStr1_CompareStringLengths__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_CompareStringLengths.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadStr2_CompareStringLengths__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_CompareStringLengths.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}

void ProgramSynthesisExperiment::Inst_LoadStr3_CompareStringLengths__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_CompareStringLengths.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[2]);
}

void ProgramSynthesisExperiment::Inst_SubmitVal_CompareStringLengths__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_CompareStringLengths.Submit((bool)wmem.AccessVal(posA).GetNum());
}

// Numeric args
void ProgramSynthesisExperiment::Inst_LoadStr1_CompareStringLengths__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_CompareStringLengths.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadStr2_CompareStringLengths__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_CompareStringLengths.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}

void ProgramSynthesisExperiment::Inst_LoadStr3_CompareStringLengths__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_CompareStringLengths.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[2]);
}

void ProgramSynthesisExperiment::Inst_SubmitVal_CompareStringLengths__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_CompareStringLengths.Submit((bool)wmem.AccessVal(posA).GetNum());
} 

// with/without type searching
void ProgramSynthesisExperiment::Inst_SubmitVal_CompareStringLengths__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_CompareStringLengths.Submit((bool)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitVal_CompareStringLengths__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_nums[0], hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_CompareStringLengths.Submit((bool)wmem.AccessVal(posA).GetNum());
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// -- String length backwards --
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tag-based args
void ProgramSynthesisExperiment::Inst_LoadStrVec_StringLengthsBackwards__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_StringLengthsBackwards.GetTestInput(eval_util.use_training_set, eval_util.current_testID));
}

void ProgramSynthesisExperiment::Inst_SubmitVal_StringLengthsBackwards__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_StringLengthsBackwards.Submit((size_t)wmem.AccessVal(posA).GetNum());
}

// Numeric args
void ProgramSynthesisExperiment::Inst_LoadStrVec_StringLengthsBackwards__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_StringLengthsBackwards.GetTestInput(eval_util.use_training_set, eval_util.current_testID));
}

void ProgramSynthesisExperiment::Inst_SubmitVal_StringLengthsBackwards__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_StringLengthsBackwards.Submit((size_t)wmem.AccessVal(posA).GetNum());
}

// with/without type searching
void ProgramSynthesisExperiment::Inst_SubmitVal_StringLengthsBackwards__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_StringLengthsBackwards.Submit((size_t)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitVal_StringLengthsBackwards__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_nums[0], hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_StringLengthsBackwards.Submit((size_t)wmem.AccessVal(posA).GetNum());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// -- Last index of zero --
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// tag-based args
void ProgramSynthesisExperiment::Inst_LoadVec_LastIndexOfZero__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_LastIndexOfZero.GetTestInput(eval_util.use_training_set, eval_util.current_testID));
}
// numeric args
void ProgramSynthesisExperiment::Inst_LoadVec_LastIndexOfZero__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_LastIndexOfZero.GetTestInput(eval_util.use_training_set, eval_util.current_testID));
}

void ProgramSynthesisExperiment::Inst_SubmitNum_LastIndexOfZero__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_LastIndexOfZero.Submit((int)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_LastIndexOfZero__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_LastIndexOfZero.Submit((int)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_LastIndexOfZero__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_LastIndexOfZero.Submit((int)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitNum_LastIndexOfZero__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_nums[0], hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_LastIndexOfZero.Submit((int)wmem.AccessVal(posA).GetNum());
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// -- Mirror image --
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// no args
void ProgramSynthesisExperiment::Inst_SubmitTrue_MirrorImage(hardware_t & hw, const inst_t & inst) { prob_utils_MirrorImage.Submit(true); }
void ProgramSynthesisExperiment::Inst_SubmitFalse_MirrorImage(hardware_t & hw, const inst_t & inst) { prob_utils_MirrorImage.Submit(false); }
// tag args
void ProgramSynthesisExperiment::Inst_LoadVec1_MirrorImage__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_MirrorImage.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadVec2_MirrorImage__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_MirrorImage.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}

// num args
void ProgramSynthesisExperiment::Inst_LoadVec1_MirrorImage__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_MirrorImage.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadVec2_MirrorImage__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_MirrorImage.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}

// with and without type searching
void ProgramSynthesisExperiment::Inst_SubmitVal_MirrorImage__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_MirrorImage.Submit((bool)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitVal_MirrorImage__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_MirrorImage.Submit((bool)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitVal_MirrorImage__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::NUM)) return;

  prob_utils_MirrorImage.Submit((bool)wmem.AccessVal(posA).GetNum());
}

void ProgramSynthesisExperiment::Inst_SubmitVal_MirrorImage__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_nums[0], hardware_t::MemPosType::NUM);
  if (!hw.IsValidMemPos(posA)) return;

  prob_utils_MirrorImage.Submit((bool)wmem.AccessVal(posA).GetNum());
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// -- Vectors summed --  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// tag args
void ProgramSynthesisExperiment::Inst_LoadVec1_VectorsSummed__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_VectorsSummed.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadVec2_VectorsSummed__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_VectorsSummed.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}

// num args
void ProgramSynthesisExperiment::Inst_LoadVec1_VectorsSummed__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_MirrorImage.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[0]);
}

void ProgramSynthesisExperiment::Inst_LoadVec2_VectorsSummed__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(posA)) return;

  wmem.Set(posA, prob_utils_MirrorImage.GetTestInput(eval_util.use_training_set, eval_util.current_testID)[1]);
}

// with and without type searching
void ProgramSynthesisExperiment::Inst_SubmitVec_VectorsSummed__TAG_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity(), hardware_t::MemPosType::VEC);
  if (!hw.IsValidMemPos(posA)) return;

  const emp::vector<hardware_t::MemoryValue> & vec = wmem.AccessVec(posA);
  emp::vector<int> output;
  for (size_t i = 0; i < vec.size(); ++i) {
    if (vec[i].GetType() == hardware_t::MemoryValue::MemoryType::NUM) {
      output.emplace_back((int)vec[i].GetNum());
    }
  }
  prob_utils_VectorsSummed.Submit(output); 
}

void ProgramSynthesisExperiment::Inst_SubmitVec_VectorsSummed__TAG_ARGS_NO_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_tags[0], hw.GetMinTagSpecificity());
  if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::VEC)) return;

  const emp::vector<hardware_t::MemoryValue> & vec = wmem.AccessVec(posA);
  emp::vector<int> output;
  for (size_t i = 0; i < vec.size(); ++i) {
    if (vec[i].GetType() == hardware_t::MemoryValue::MemoryType::NUM) {
      output.emplace_back((int)vec[i].GetNum());
    }
  }
  prob_utils_VectorsSummed.Submit(output); 
}

void ProgramSynthesisExperiment::Inst_SubmitVec_VectorsSummed__NUM_ARGS(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  emp_assert(inst.arg_nums.size() >= 1);
  size_t posA = inst.arg_nums[0]; if (!hw.IsValidMemPos(wmem, posA, hardware_t::MemPosType::VEC)) return;

  const emp::vector<hardware_t::MemoryValue> & vec = wmem.AccessVec(posA);
  emp::vector<int> output;
  for (size_t i = 0; i < vec.size(); ++i) {
    if (vec[i].GetType() == hardware_t::MemoryValue::MemoryType::NUM) {
      output.emplace_back((int)vec[i].GetNum());
    }
  }
  prob_utils_VectorsSummed.Submit(output); 
}

void ProgramSynthesisExperiment::Inst_SubmitVec_VectorsSummed__NUM_ARGS_WITH_TYPE_SEARCH(hardware_t & hw, const inst_t & inst) {
  hardware_t::CallState & state = hw.GetCurCallState();
  hardware_t::Memory & wmem = state.GetWorkingMem();

  // Find arguments
  size_t posA = hw.FindBestMemoryMatch(wmem, inst.arg_nums[0], hardware_t::MemPosType::VEC);
  if (!hw.IsValidMemPos(posA)) return;

  const emp::vector<hardware_t::MemoryValue> & vec = wmem.AccessVec(posA);
  emp::vector<int> output;
  for (size_t i = 0; i < vec.size(); ++i) {
    if (vec[i].GetType() == hardware_t::MemoryValue::MemoryType::NUM) {
      output.emplace_back((int)vec[i].GetNum());
    }
  }
  prob_utils_VectorsSummed.Submit(output); 
}

#endif