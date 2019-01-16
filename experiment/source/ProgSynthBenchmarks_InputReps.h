#ifndef PROG_SYNTH_BENCHMARK_INPUT_REPRESENTATIONS_H
#define PROG_SYNTH_BENCHMARK_INPUT_REPRESENTATIONS_H

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <unordered_set>

#include "base/Ptr.h"
#include "base/vector.h"
#include "tools/Random.h"
#include "tools/random_utils.h"
#include "tools/math.h"
#include "tools/vector_utils.h"
#include "tools/sequence_utils.h"
#include "tools/string_utils.h"
#include "tools/stats.h"

#include "parser.hpp"

#include "TestCaseSet.h"
#include "Utilities.h"

/*

Genetic representations for different input types used across the programming
synthesis benchmarks.

Copied over from antagonistic lexicase selection work.

*/

// Problem-specific, static variables

// -- PROBLEM: Small or Large --
constexpr int SmallOrLarge__SMALL_THRESH = 1000;
constexpr int SmallOrLarge__LARGE_THRESH = 2000;
const std::string SmallOrLarge__SMALL_STR = "small";
const std::string SmallOrLarge__LARGE_STR = "large";
const std::string SmallOrLarge__NONE_STR = "";

const std::string Grade__A_STR = "A";
const std::string Grade__B_STR = "B";
const std::string Grade__C_STR = "C";
const std::string Grade__D_STR = "D";
const std::string Grade__F_STR = "F";

// ================== Problem I/O Type Aliases ==================
using Problem_NumberIO_input_t = std::pair<int, double>;
using Problem_NumberIO_output_t = double;

using Problem_SmallOrLarge_input_t = int;
using Problem_SmallOrLarge_output_t = std::string;

using Problem_ForLoopIndex_input_t = std::array<int, 3>;
using Problem_ForLoopIndex_output_t = emp::vector<int>;

using Problem_CompareStringLengths_input_t = std::array<std::string, 3>;
using Problem_CompareStringLengths_output_t = bool;

using Problem_DoubleLetters_input_t = std::string;
using Problem_DoubleLetters_output_t = std::string;

using Problem_CollatzNumbers_input_t = int;
using Problem_CollatzNumbers_output_t = int;

using Problem_ReplaceSpaceWithNewline_input_t = std::string;
using Problem_ReplaceSpaceWithNewline_output_t = std::pair<std::string, int>;

using Problem_StringDifferences_input_t = std::array<std::string, 2>;
using Problem_StringDifferences_output_t = std::string;

using Problem_EvenSquares_input_t = int;
using Problem_EvenSquares_output_t = int;

using Problem_WallisPi_input_t = int;
using Problem_WallisPi_output_t = double;

using Problem_StringLengthsBackwards_input_t = emp::vector<std::string>;
using Problem_StringLengthsBackwards_output_t = emp::vector<size_t>;

using Problem_LastIndexOfZero_input_t = emp::vector<int>;
using Problem_LastIndexOfZero_output_t = int;

using Problem_VectorAverage_input_t = emp::vector<double>;
using Problem_VectorAverage_output_t = double;

using Problem_CountOdds_input_t = emp::vector<int>;
using Problem_CountOdds_output_t = int;

using Problem_MirrorImage_input_t = std::array<emp::vector<int>, 2>;
using Problem_MirrorImage_output_t = bool;

using Problem_SuperAnagrams_input_t = std::array<std::string, 2>;
using Problem_SuperAnagrams_output_t = bool;

using Problem_SumOfSquares_input_t = int;
using Problem_SumOfSquares_output_t = int;

using Problem_VectorsSummed_input_t = std::array<emp::vector<int>, 2>;
using Problem_VectorsSummed_output_t = emp::vector<int>;

using Problem_XWordLines_input_t = std::pair<int, std::string>;
using Problem_XWordLines_output_t = std::string;

using Problem_PigLatin_input_t = std::string;
using Problem_PigLatin_output_t = std::string;

using Problem_NegativeToZero_input_t = emp::vector<int>;
using Problem_NegativeToZero_output_t = emp::vector<int>;

using Problem_ScrabbleScore_input_t = std::string;
using Problem_ScrabbleScore_output_t = int;

using Problem_Checksum_input_t = std::string;
using Problem_Checksum_output_t = std::string;

using Problem_Digits_input_t = int;
using Problem_Digits_output_t = emp::vector<int>;

using Problem_Grade_input_t = std::array<int, 5>;
using Problem_Grade_output_t = std::string;

using Problem_Median_input_t = std::array<int, 3>;
using Problem_Median_output_t = int;

using Problem_Smallest_input_t = std::array<int, 4>;
using Problem_Smallest_output_t = int;

using Problem_Syllables_input_t = std::string;
using Problem_Syllables_output_t = std::string;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: NumberIO
// - Input type: [double, integer]
// - Output type: double 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Problem_NumberIO_input_t GenRandomTestInput_NumberIO(emp::Random & rand, const std::pair<int, int> & int_range, const std::pair<double, double> & double_range) {
  emp_assert(double_range.first < double_range.second);
  emp_assert(int_range.first < int_range.second);
  return Problem_NumberIO_input_t{rand.GetInt(int_range.first, int_range.second), rand.GetDouble(double_range.first, double_range.second)};
}
  
Problem_NumberIO_output_t GenCorrectOut_NumberIO(const Problem_NumberIO_input_t & input) {
  return input.first + input.second;
}

struct ProblemUtilities_NumberIO {
  using input_t = Problem_NumberIO_input_t;
  using output_t = Problem_NumberIO_output_t;
  
  using testcase_set_t = TestCaseSet<input_t,output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  // A few useful things for use within a test evaluation
  bool submitted;
  double submitted_val; // if going to do string thing, we can have a submission_str.
  double MAX_ERROR;

  ProblemUtilities_NumberIO() 
    : testing_set(ProblemUtilities_NumberIO::LoadTestCaseFromLine),
      training_set(ProblemUtilities_NumberIO::LoadTestCaseFromLine),
      submitted(false), submitted_val(0.0), MAX_ERROR(1)
  { ; }

  ~ProblemUtilities_NumberIO() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }
  
  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }

  void ResetTestEval() {
    submitted = false;
    submitted_val = 0.0;
  }

  void Submit(double val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const std::string & line) {
    emp::vector<std::string> split_line = emp::slice(line, ',');
    input_t input;
    output_t output;
    input.second = std::atof(split_line[0].c_str());
    input.first = std::atof(split_line[1].c_str());
    output = std::atof(split_line[2].c_str());
    return {input, output};
  }

  /// Calculate pass/fail score on NumberIO problem.
  std::pair<double, bool> CalcScorePassFail(const Problem_NumberIO_output_t & correct_test_output, double sub) {
    const bool pass = (double)sub == correct_test_output;
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const Problem_NumberIO_output_t & correct_test_output, double sub) {
    emp_assert(MAX_ERROR != 0);
    // If output is correct, return a score of 1.0 and mark that submission passes.
    if (correct_test_output == sub) {
      return {1.0, true}; 
    } else { // Otherwise, return {score=[0:1], false}
      double error = emp::Abs(correct_test_output - sub);
      emp_assert(error != 0, "Error shouldn't equal zero here");
      double score = (error <= MAX_ERROR) ? 1 - (error/MAX_ERROR) : 0;
      return {score, false};
    }  
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << in.first << "," << in.second;
    os << "\"";
  }

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: Small or large
// - Input: integer
// - Output: string {'small', 'large', ''}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Generate random test input
Problem_SmallOrLarge_input_t GenRandomTestInput_SmallOrLarge(emp::Random & rand, const std::pair<int,int> & int_range) {
  emp_assert(int_range.first < int_range.second);
  return rand.GetInt(int_range.first, int_range.second+1);
}

Problem_SmallOrLarge_output_t GenCorrectOut_SmallOrLarge(const Problem_SmallOrLarge_input_t & input) {
  if (input < SmallOrLarge__SMALL_THRESH) return SmallOrLarge__SMALL_STR;
  else if (input >= SmallOrLarge__LARGE_THRESH) return SmallOrLarge__LARGE_STR;
  else return SmallOrLarge__NONE_STR;
}

struct ProblemUtilities_SmallOrLarge { 
  using this_t = ProblemUtilities_SmallOrLarge;
  using input_t = Problem_SmallOrLarge_input_t;
  using output_t = Problem_SmallOrLarge_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  bool submitted;
  std::string submitted_val;

  ProblemUtilities_SmallOrLarge()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val("")
  { ; }

  ~ProblemUtilities_SmallOrLarge() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }

  void ResetTestEval() {
    submitted = false;
    submitted_val = "";
  }

  void Submit(const std::string & val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const std::string & line) {
    emp::vector<std::string> split_line = emp::slice(line + " ", ',');
    input_t input;    // int
    output_t output;  // std::string
    // std::cout << "LINE=" << line << std::endl;
    input = std::atof(split_line[0].c_str());
    output = split_line[1];
    if (!(output == " " || output == "small " || output == "large ")) {
      std::cout << "ERROR! Bad output ("<<output<<") from line: " << line << std::endl;
      exit(-1);
    }
    if (output == " ") output = "";
    else if (output == "small ") output = "small";
    else if (output == "large ") output = "large";
    return {input, output};
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << in;
    os << "\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: ForLoopIndex
// - Input: std::array<int, 3>;
// - Output: emp::vector<int>;
// - Description: Given 3 integer inputs (start, end, step), print the integers in the sequence:
//              - n0 = start
//              - ni = ni-1 + step
//              - for each ni < end
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_ForLoopIndex_input_t GenRandomTestInput_ForLoopIndex(emp::Random & rand,
                                                         const std::pair<int,int> & start_end_range,
                                                         const std::pair<int,int> & step_range) {
  // Guarantee that start comes before end.
  int start, end, step;
  start = rand.GetInt(start_end_range.first, start_end_range.second+1);
  end = rand.GetInt(start_end_range.first, start_end_range.second+1);
  step = rand.GetInt(step_range.first, step_range.second+1);

  // Need to follow the following rules:
  // (1) start < end               // -> Start should come before end
  // (2) start + (20xstep)+1 > end // -> Enumeration should take no more than 20 steps
  while (true) {
    if (end < start) std::swap(end, start);
    if (start + (20*step) + 1 > end) break;
    start = rand.GetInt(start_end_range.first, start_end_range.second+1);
    end = rand.GetInt(start_end_range.first, start_end_range.second+1);
    step = rand.GetInt(step_range.first, step_range.second+1);
  }
  return {start, end, step};
}

/// Generate correct output
Problem_ForLoopIndex_output_t GenCorrectOut_ForLoopIndex(const Problem_ForLoopIndex_input_t & input) {
  emp::vector<int> out;
  for (int i = input[0]; i < input[1]; i+= input[2]) {
    out.emplace_back(i);
  }
  return out;
}

struct ProblemUtilities_ForLoopIndex { 
  using this_t = ProblemUtilities_ForLoopIndex;
  using input_t = Problem_ForLoopIndex_input_t;
  using output_t = Problem_ForLoopIndex_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  // --- Useful during a test evaluation ---
  bool submitted;
  emp::vector<int> submitted_val;

  ProblemUtilities_ForLoopIndex()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val()
  { ; }

  ~ProblemUtilities_ForLoopIndex() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }

  void ResetTestEval() {
    submitted = false;
    submitted_val.clear();
  }

  void Submit(int val) {
    submitted = true;
    submitted_val.emplace_back(val);
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const std::string & line) {
    emp::vector<std::string> split_line = emp::slice(line, ',');
    input_t input;    // int
    output_t output;  // std::string
    // Start = line[0]
    input[0] = std::atof(split_line[0].c_str());
    // End = line[1]
    input[1] = std::atof(split_line[1].c_str());
    // Step = line[2]
    input[2] = std::atof(split_line[2].c_str());
    output = GenCorrectOut_ForLoopIndex(input);  
    return {input, output};
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    const double max_dist = emp::Max(correct_test_output.size(), sub.size());
    double dist = emp::calc_edit_distance(correct_test_output, sub);
    if (dist == 0) {
      return {1.0, true};
    } else {
      return {(max_dist - dist)/max_dist, false};
    }
  } // todo - test this

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << in[0] << "," << in[1] << "," << in[2];
    os << "\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: CompareStringLengths
// - Input: std::array<std::string, 3>;
// - Output: bool;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// /// Generate random test input
Problem_CompareStringLengths_input_t GenRandomTestInput_CompareStringLengths(emp::Random & rand,
                                                                             std::pair<size_t, size_t> str_size_range) {
  // Valid characters: \n, \t, [32, 127)
  emp::vector<char> valid_chars = {'\n', '\t'};
  for (size_t i = 32; i < 127; ++i) valid_chars.emplace_back((char)i);

  // String 1
  size_t str_size = rand.GetUInt(str_size_range.first, str_size_range.second+1);
  std::string str1;
  for (size_t i = 0; i < str_size; ++i) str1.push_back(valid_chars[rand.GetUInt(valid_chars.size())]);

  // String 2
  str_size = rand.GetUInt(str_size_range.first, str_size_range.second+1);
  std::string str2;
  for (size_t i = 0; i < str_size; ++i) str2.push_back(valid_chars[rand.GetUInt(valid_chars.size())]);

  // String 3
  str_size = rand.GetUInt(str_size_range.first, str_size_range.second+1);
  std::string str3;
  for (size_t i = 0; i < str_size; ++i) str3.push_back(valid_chars[rand.GetUInt(valid_chars.size())]);
  // std::cout << "---- RANDOM STRING ----" << std::endl;
  // std::cout << "String 1: " << str1 << std::endl;
  // std::cout << "String 2: " << str2 << std::endl;
  // std::cout << "String 3: " << str3 << std::endl;
  return {str1, str2, str3};
}

/// Generate correct output
Problem_CompareStringLengths_output_t GenCorrectOut_CompareStringLengths(const Problem_CompareStringLengths_input_t & input) {
  if (input[0].size() < input[1].size() && input[1].size() < input[2].size()) return true;
  return false;
}

struct ProblemUtilities_CompareStringLengths { 
  using this_t = ProblemUtilities_CompareStringLengths;
  using input_t = Problem_CompareStringLengths_input_t;
  using output_t = Problem_CompareStringLengths_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  // --- Useful during a test evaluation ---
  bool submitted;
  bool submitted_val;

  ProblemUtilities_CompareStringLengths()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(false)
  { ; }

  ~ProblemUtilities_CompareStringLengths() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }

  void ResetTestEval() {
    submitted = false;
    submitted_val = false;
  }

  void Submit(bool val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    // Load input.
    input[0] = line[0];
    input[1] = line[1];
    input[2] = line[2];
    // Load output.
    if (line[3] == "false") output = false;
    else if (line[3] == "true") output = true;
    else {
      std::cout << "ERROR when loading test case output (" << line[3] << "). Exiting." << std::endl;
      exit(-1);
    }
    return {input, output};
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"[{BEGIN-STR}" << in[0] << "{END-STR},{BEGIN-STR}" << in[1] << "{END-STR},{BEGIN-STR}" << in[2] << "{END-STR}]\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: LastIndexOfZero
// - Input: emp::vector<int>; (at least one value *must* be zero)
// - Output: int;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_LastIndexOfZero_input_t GenRandomTestInput_LastIndexOfZero(emp::Random & rand,
                                                                          const std::pair<size_t, size_t> & vec_size_range,
                                                                          const std::pair<int, int> & vec_val_range) {
  emp::vector<int> input;
  size_t zero_cnt = 0;
  // Randomize size.
  const size_t vec_size = rand.GetUInt(vec_size_range.first, vec_size_range.second+1);
  for (size_t i = 0; i < vec_size; ++i) {
    int val = rand.GetInt(vec_val_range.first, vec_val_range.second+1);
    input.emplace_back(val);
    if (val == 0) zero_cnt++;
  }
  // ensure there's at least one zero
  if (zero_cnt == 0) {
    input[rand.GetUInt(0, input.size())] = 0;
  }
  emp_assert(emp::Has(input, 0));
  return input;
}

Problem_LastIndexOfZero_output_t GenCorrectOut_LastIndexOfZero(const Problem_LastIndexOfZero_input_t & input) {
  emp_assert(emp::Has(input, 0));
  for (int i = input.size()-1; i >= 0; --i) {
    if (input[i] == 0) return i;
  }
  std::cout << "ERROR! Failed to find 0 in last index of zero output!" << std::endl;
  exit(-1);
  return 0; // Should never get here!
}

struct ProblemUtilities_LastIndexOfZero { 
  using this_t = ProblemUtilities_LastIndexOfZero;
  using input_t = Problem_LastIndexOfZero_input_t;
  using output_t = Problem_LastIndexOfZero_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  // --- Useful during a test evaluation ---
  bool submitted;
  int submitted_val;

  int MAX_ERROR;
  
  ProblemUtilities_LastIndexOfZero()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0)
  { ; }

  ~ProblemUtilities_LastIndexOfZero() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }


  void ResetTestEval() {
    submitted = false;
    submitted_val = 0;
  }

  void Submit(int val) {
    submitted = true;
    submitted_val = val;
  }


  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    
    std::string input_str = line[0];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }
    emp::vector<std::string> sliced_input_str = emp::slice(input_str, ' ');
    
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input.emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }
    
    // Calculate correct output given loaded input
    int calc_out = GenCorrectOut_LastIndexOfZero(input);
    // Get output from file
    output = std::atoi(line[1].c_str());
    // make sure generated output and read output match
    if (calc_out != output) {
      std::cout << "ERROR! Generated last index of zero output does not match read output! Exiting." << std::endl;
      exit(-1);
    }

    return {input, output};
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    if (correct_test_output == sub) {
      return {1.0, true};
    } else {
      double error = (double)emp::Abs(correct_test_output - sub);
      emp_assert(error != 0, "Error shouldn't be zero here!");
      double score = (error <= MAX_ERROR) ? 1 - (error/(double)MAX_ERROR) : 0.0;
      return {score, false};
    }
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << "[";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << in[i];
    }
    os << "]"; 
    os << "\"";
  }

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: MirrorImage
// - Input: std::array<emp::vector<int>, 2>;
// - Output: bool;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random input
/// Generate random test input
Problem_MirrorImage_input_t GenRandomTestInput_MirrorImage(emp::Random & rand,
                                                           const std::pair<size_t, size_t> & vec_size_range,
                                                           const std::pair<int, int> & vec_val_range) {
  std::array<emp::vector<int>, 2> input;
  // 4 random cases: equal, random, mirrored, nearly mirrored
  const size_t vec_size = rand.GetUInt(vec_size_range.first, vec_size_range.second+1);
  const size_t rand_case = rand.GetUInt(0, 4);
  switch (rand_case) {
    case 0: { // Equal
      for (size_t k = 0; k < vec_size; ++k) {
        input[0].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
      }
      input[1] = input[0];
      break; 
    }
    case 1: { // Random
      for (size_t i = 0; i < input.size(); ++i) {
        for (size_t k = 0; k < vec_size; ++k) {
          input[i].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
        }
      } 
      break; 
    }
    case 2: { // Mirrored
      for (size_t k = 0; k < vec_size; ++k) {
        input[0].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
      }
      input[1] = input[0];
      std::reverse(std::begin(input[1]), std::end(input[1]));
      break; 
    }
    case 3: { // Nearly mirrored
      for (size_t k = 0; k < vec_size; ++k) {
        input[0].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
      }
      input[1] = input[0];
      std::reverse(std::begin(input[1]), std::end(input[1]));
      const size_t num_randos = rand.GetUInt(input[0].size());
      for (size_t i = 0; i < num_randos; ++i) {
        input[0][rand.GetUInt(input[0].size())] = rand.GetInt(vec_val_range.first, vec_val_range.second+1);
      }
      break; 
    }
  }
  return input;
}

/// Generate correct output given input
Problem_MirrorImage_output_t GenCorrectOut_MirrorImage(const Problem_MirrorImage_input_t & input) {
  emp_assert(input[0].size() == input[1].size());
  if (input[0].size() != input[1].size()) return false;
  const size_t vec_size = input[0].size();
  for (size_t i = 0; i < vec_size; ++i) {
    const size_t ri = (vec_size - 1) - i;
    if (input[0][i] != input[1][ri]) return false;
  }
  return true;
}

struct ProblemUtilities_MirrorImage {
  using this_t = ProblemUtilities_MirrorImage;
  using input_t = Problem_MirrorImage_input_t;
  using output_t = Problem_MirrorImage_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  // --- Useful during a test evaluation ---
  bool submitted;
  bool submitted_val;
  
  ProblemUtilities_MirrorImage()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0)
  { ; }

  ~ProblemUtilities_MirrorImage() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }


  void ResetTestEval() {
    submitted = false;
    submitted_val = false;
  }

  void Submit(bool val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    
    // Vector 1
    std::string input_str = line[0];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }
    emp::vector<std::string> sliced_input_str = emp::slice(input_str, ' ');
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input[0].emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }

    // Vector 2
    input_str = line[1];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }  
    sliced_input_str = emp::slice(input_str, ' ');
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input[1].emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }   

    // Calculate correct output given loaded input
    bool calc_out = GenCorrectOut_MirrorImage(input);
    // Get output from file
    if (line[2] == "false") output = false;
    else if (line[2] == "true") output = true;
    else {
      std::cout << "Unrecognized output value for mirror image examples! Exiting." << std::endl;
      exit(-1);
    }

    // make sure generated output and read output match
    if (calc_out != output) {
      std::cout << "ERROR! Generated output does not match read output! Exiting." << std::endl;
      exit(-1);
    }

    return {input, output};
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << "[[";
    for (size_t i = 0; i < in[0].size(); ++i) {
      if (i) os << ",";
      os << in[0][i];
    }
    os << "],["; 
    for (size_t i = 0; i < in[1].size(); ++i) {
      if (i) os << ",";
      os << in[1][i];
    }
    os << "]]"; 
    os << "\"";
  }

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: SumOfSquares (NOTE: only for: 1 >= n <= 100)
// - Input: int
// - Output: int
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::array<int, 101> SumOfSquaresLookup = {0, 1, 5, 14, 30, 55, 91, 140, 204, 285, 385, 506, 650, 819, 1015, 1240, 1496, 1785, 2109, 2470, 2870, 3311, 3795, 4324, 4900, 5525, 6201, 6930, 7714, 8555, 9455, 10416, 11440, 12529, 13685, 14910, 16206, 17575, 19019, 20540, 22140, 23821, 25585, 27434, 29370, 31395, 33511, 35720, 38024, 40425, 42925, 45526, 48230, 51039, 53955, 56980, 60116, 63365, 66729, 70210, 73810, 77531, 81375, 85344, 89440, 93665, 98021, 102510, 107134, 111895, 116795, 121836, 127020, 132349, 137825, 143450, 149226, 155155, 161239, 167480, 173880, 180441, 187165, 194054, 201110, 208335, 215731, 223300, 231044, 238965, 247065, 255346, 263810, 272459, 281295, 290320, 299536, 308945, 318549, 328350, 338350};

Problem_SumOfSquares_input_t GenRandomTestInput_SumOfSquares(emp::Random & rnd, const std::pair<int, int> & num_range) {
  return rnd.GetInt(num_range.first, num_range.second+1);
}

Problem_SumOfSquares_output_t GenCorrectOut_SumOfSquares(const Problem_SumOfSquares_input_t & input) {
  if (input < (int)SumOfSquaresLookup.size()) {
    return SumOfSquaresLookup[(size_t)input];
  } else {
    return ((input) * (input + 1) * ((2*input)+1))/6;
  }
}

struct ProblemUtilities_SumOfSquares { 
  using this_t = ProblemUtilities_SumOfSquares;
  using input_t = Problem_SumOfSquares_input_t;
  using output_t = Problem_SumOfSquares_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  // // --- Useful during a test evaluation ---
  bool submitted;
  int submitted_val;
  int MAX_ERROR;

  ProblemUtilities_SumOfSquares()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0)
  { ; }

  ~ProblemUtilities_SumOfSquares() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }


  void ResetTestEval() {
    submitted = false;
    submitted_val = 0;
  }

  void Submit(int val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    // Load input.
    input = std::atof(line[0].c_str());
    // Load output.
    output = std::atof(line[1].c_str());
    emp_assert(output == GenCorrectOut_SumOfSquares(input));
    return {input, output};
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    if (correct_test_output == sub) {
      return {1.0, true};
    } else {
      double error = (double)emp::Abs(correct_test_output - sub);
      emp_assert(error != 0, "Error shouldn't be zero here!");
      double score = (error <= MAX_ERROR) ? 1 - (error/(double)MAX_ERROR) : 0.0;
      return {score, false};
    }
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << in;
    os << "\"";
  }

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: VectorsSummed
// - Input: std::array<emp::vector<int>, 2>;
// - Output: emp::vector<int>;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_VectorsSummed_input_t GenRandomTestInput_VectorsSummed(emp::Random & rand,
                                                               const std::pair<size_t, size_t> & vec_size_range,
                                                                const std::pair<int, int> & vec_val_range) {
  std::array<emp::vector<int>, 2> input; // Input
  const size_t vec_size = rand.GetUInt(vec_size_range.first, vec_size_range.second+1);
  for (size_t i = 0; i < vec_size; ++i) {
    input[0].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
    input[1].emplace_back(rand.GetInt(vec_val_range.first, vec_val_range.second+1));
  }
  return input;
}

/// Generate correct output given input
Problem_VectorsSummed_output_t GenCorrectOut_VectorsSummed(const Problem_VectorsSummed_input_t & input) {
  emp::vector<int> output;
  emp_assert(input[0].size() == input[1].size());
  const size_t vec_size = input[0].size();
  for (size_t i = 0; i < vec_size; ++i) {
    const int res = input[0][i] + input[1][i];
    output.emplace_back(res);
  }
  emp_assert(output.size() == input[0].size());
  return output;
}

struct ProblemUtilities_VectorsSummed { 
  using this_t = ProblemUtilities_VectorsSummed;
  using input_t = Problem_VectorsSummed_input_t;
  using output_t = Problem_VectorsSummed_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  // --- Useful during a test evaluation ---
  bool submitted;
  emp::vector<int> submitted_val;
  int MAX_NUM;
  int MAX_ERROR;
  
  ProblemUtilities_VectorsSummed()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0)
  { ; }

  ~ProblemUtilities_VectorsSummed() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }


  void ResetTestEval() {
    submitted = false;
    submitted_val.clear();
  }

  void Submit(const emp::vector<int> & vec) {
    submitted = true;
    submitted_val = vec;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    
    // Vector 1
    std::string input_str = line[0];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }
    emp::vector<std::string> sliced_input_str = emp::slice(input_str, ' ');
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input[0].emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }

    // Vector 2
    input_str = line[1];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }  
    sliced_input_str = emp::slice(input_str, ' ');
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      input[1].emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }   

    // Output vector
    input_str = line[2];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }  
    sliced_input_str = emp::slice(input_str, ' ');
    for (size_t i = 0; i < sliced_input_str.size(); ++i) {
      output.emplace_back(std::atoi(sliced_input_str[i].c_str()));
    }

    // Calculate correct output given loaded input
    emp::vector<int> calc_out = GenCorrectOut_VectorsSummed(input);
    
    // make sure generated output and read output match
    if (calc_out != output) {
      std::cout << "ERROR! Generated output does not match read output! Exiting." << std::endl;
      exit(-1);
    }

    return {input, output};
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    if (correct_test_output == sub) {
      return {1.0, true};
    } else {
      // double error = (double)emp::Abs(correct_test_output - sub);
      double error = 0;
      for (size_t i = 0; i < correct_test_output.size(); ++i) {
        if (i < sub.size()) {
          // Add error.
          error += emp::Abs((double)((double)correct_test_output[i] - (double)sub[i]));
        } else {
          // Add max error.
          error += (2*MAX_NUM);
        }
      }
      // Add error for each extra thing in sub
      if (sub.size() > correct_test_output.size()) {
        error += (2*MAX_NUM) * (sub.size() - correct_test_output.size());
      }
      // Make sure programs that overflow error to be negative (WTF?) don't get really high scores.
      double score = (error <= MAX_ERROR && error >= 0) ? 1 - (error/(double)MAX_ERROR) : 0.0;

      return {score, false};
    }
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << "[[";
    for (size_t i = 0; i < in[0].size(); ++i) {
      if (i) os << ",";
      os << in[0][i];
    }
    os << "],["; 
    for (size_t i = 0; i < in[1].size(); ++i) {
      if (i) os << ",";
      os << in[1][i];
    }
    os << "]]"; 
    os << "\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: Grade
// - Input: Array<Integer, 5>
// - Output: String
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Generate random test input
Problem_Grade_input_t GenRandomTestInput_Grade(emp::Random & rand, const std::pair<int, int> & num_range) {
  emp_assert(num_range.first < num_range.second);
  emp_assert(num_range.first > -1);
  Problem_Grade_input_t input;
  input[0] = -1;
  input[1] = -1;
  input[2] = -1;
  input[3] = -1;
  input[4] = -1;

  for (size_t i = 0; i < 4; ++i) {
    int val = rand.GetInt(num_range.first+1, num_range.second);
    // while (emp::Has(input, val)) val = rand.GetInt(num_range.first+1, num_range.second);
    while (input[0] == val || input[1] == val || input[2] == val || input[3] == val) val = rand.GetInt(num_range.first+1, num_range.second);
    input[i] = val;
  }
  // Sort input
  std::sort(input.begin(), input.end());
  std::reverse(std::begin(input), std::end(input));
  input[4] = rand.GetInt(num_range.first, num_range.second+1); // Grade?

  emp_assert(100 >= input[0], input[0], num_range.second);
  emp_assert(num_range.second >= input[0]);
  emp_assert(input[0] > input[1]);
  emp_assert(input[1] > input[2]);
  emp_assert(input[2] > input[3]);
  emp_assert(input[3] >= num_range.first);
  emp_assert(input[3] >= 0);
  return input;
}

// Generate correct output for a given test
Problem_Grade_output_t GenCorrectOut_Grade(const Problem_Grade_input_t & input) {
  if (input[4] >= input[0]) { return Grade__A_STR; }
  else if (input[4] >= input[1]) { return Grade__B_STR; }
  else if (input[4] >= input[2]) { return Grade__C_STR; }
  else if (input[4] >= input[3]) { return Grade__D_STR; }
  else { return Grade__F_STR; }
}

struct ProblemUtilities_Grade { 
  using this_t = ProblemUtilities_Grade;
  using input_t = Problem_Grade_input_t;
  using output_t = Problem_Grade_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  // // --- Useful during a test evaluation ---
  bool submitted;
  std::string submitted_val;

  ProblemUtilities_Grade()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val("")
  { ; }

  ~ProblemUtilities_Grade() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }


  void ResetTestEval() {
    submitted = false;
    submitted_val = "";
  }

  void Submit(const std::string & val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;
    output_t output;

    // Load input
    input[0] = std::atof(line[0].c_str());
    input[1] = std::atof(line[1].c_str());
    input[2] = std::atof(line[2].c_str());
    input[3] = std::atof(line[3].c_str());
    input[4] = std::atof(line[4].c_str());
    
    // std::cout << "=== Loading test case from file ===" << std::endl;
    // std::cout << "A thresh: " << input[0] << std::endl;
    // std::cout << "B thresh: " << input[1] << std::endl;
    // std::cout << "C thresh: " << input[2] << std::endl;
    // std::cout << "D thresh: " << input[3] << std::endl;
    // std::cout << "Grade: " << input[4] << std::endl;

    if (line[5] == "Student has a A grade.") {
      output = Grade__A_STR;
    } else if (line[5] == "Student has a B grade.") {
      output = Grade__B_STR;
    } else if (line[5] == "Student has a C grade.") {
      output = Grade__C_STR;
    } else if (line[5] == "Student has a D grade.") {
      output = Grade__D_STR;
    } else if (line[5] == "Student has a F grade.") {
      output = Grade__F_STR;
    } else {
      std::cout << "ERROR ERROR! OH NO! INVALID OUTPUT FROM GRADE EXAMPLES!" << std::endl;
      exit(-1);
    }

    output_t gen_out = GenCorrectOut_Grade(input);
    // std::cout << "Output: " << gen_out << ", " << output << std::endl;
    emp_assert(gen_out == output);

    emp_assert(100 >= input[0]);
    emp_assert(input[0] > input[1]);
    emp_assert(input[1] > input[2]);
    emp_assert(input[2] > input[3]);
    emp_assert(input[3] >= 0);

    return {input, output};

  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << in[i];
    }
    os << "\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: Median
// - Input: std::array<int, 3>;
// - Output: int;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_Median_input_t GenRandomTestInput_Median(emp::Random & rand,
                                                 const std::pair<int, int> & num_range) {
  Problem_Median_input_t input;
  for (size_t i = 0; i < input.size(); ++i) {
    input[i] = rand.GetInt(num_range.first, num_range.second+1);
  }
  return input;
}

/// Generate correct output
Problem_Median_output_t GenCorrectOut_Median(const Problem_Median_input_t & input) {
  const int min_val = emp::Min(input[0], input[1], input[2]);
  const int max_val = emp::Max(input[0], input[1], input[2]);
  // std::cout << "Input: " << input[0] << ", " << input[1] << ", " << input[2] << "; Output: " << (input[0] + input[1] + input[2]) - min_val - max_val << std::endl;
  return (input[0] + input[1] + input[2]) - min_val - max_val;
}

struct ProblemUtilities_Median { 
  using this_t = ProblemUtilities_Median;
  using input_t = Problem_Median_input_t;
  using output_t = Problem_Median_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  // // --- Useful during a test evaluation ---
  bool submitted;
  int submitted_val;

  ProblemUtilities_Median()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(0)
  { ; }

  ~ProblemUtilities_Median() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }


  void ResetTestEval() {
    submitted = false;
    submitted_val = 0;
  }

  void Submit(int val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    // Load input.
    input[0] = std::atof(line[0].c_str());
    input[1] = std::atof(line[1].c_str());
    input[2] = std::atof(line[2].c_str());
    // Load output.
    output = std::atof(line[3].c_str());
    emp_assert(output == GenCorrectOut_Median(input));
    return {input, output};
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << in[i];
    }
    os << "\"";
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: Smallest
// - Input: Array<Integer, 4>
// - Output: int
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generate random test input
Problem_Smallest_input_t GenRandomTestInput_Smallest(emp::Random & rand,
                                                     const std::pair<int, int> & num_range) {
  Problem_Smallest_input_t input;
  for (size_t i = 0; i < input.size(); ++i) {
    input[i] = rand.GetInt(num_range.first, num_range.second+1);
  }
  return input;
}

/// Generate correct output
Problem_Smallest_output_t GenCorrectOut_Smallest(const Problem_Smallest_input_t & input) {
  int smallest = input[0];
  for (size_t i = 1; i < input.size(); ++i) {
    if (input[i] < smallest) smallest = input[i];
  }
  return smallest;
}

struct ProblemUtilities_Smallest { 
  using this_t = ProblemUtilities_Smallest;
  using input_t = Problem_Smallest_input_t;
  using output_t = Problem_Smallest_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  // // --- Useful during a test evaluation ---
  bool submitted;
  int submitted_val;

  ProblemUtilities_Smallest()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val(false)
  { ; }

  ~ProblemUtilities_Smallest() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }

  void ResetTestEval() {
    submitted = false;
    submitted_val = 0;
  }

  void Submit(int val) {
    submitted = true;
    submitted_val = val;
  }

  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    // Load input.
    input[0] = std::atof(line[0].c_str());
    input[1] = std::atof(line[1].c_str());
    input[2] = std::atof(line[2].c_str());
    input[3] = std::atof(line[3].c_str());
    // Load output.
    output = std::atof(line[4].c_str());
    emp_assert(output == GenCorrectOut_Smallest(input));
    return {input, output};
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << in[i];
    }
    os << "\"";
  }

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problem: StringLengthsBackwards
// - Input: emp::vector<std::string> 
// - Output: emp::vector<size_t>
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// Generate random test input
Problem_StringLengthsBackwards_input_t GenRandomTestInput_StringLengthsBackwards(emp::Random & rand,
                                                                                 const std::pair<size_t, size_t> & str_cnt_range,
                                                                                 const std::pair<size_t, size_t> & str_size_range) {

  emp::vector<char> valid_chars = {'\n', '\t'};
  for (size_t i = 32; i < 127; ++i) valid_chars.emplace_back((char)i);

  emp::vector<std::string> rand_input;
  // How many strings should we generate?
  const size_t str_cnt = rand.GetUInt(str_cnt_range.first, str_cnt_range.second);
  // Generate each string randomly.
  for (size_t i = 0; i < str_cnt; ++i) {
    const size_t str_size = rand.GetUInt(str_size_range.first, str_size_range.second);
    rand_input.emplace_back("");
    for (size_t i = 0; i < str_size; ++i) rand_input.back().push_back(valid_chars[rand.GetUInt(valid_chars.size())]);
  }
  return rand_input;
}

/// Generate correct output
Problem_StringLengthsBackwards_output_t GenCorrectOut_StringLengthsBackwards(const Problem_StringLengthsBackwards_input_t & input) {
  emp::vector<size_t> output;
  for (int i = input.size()-1; i >= 0; --i) {
    output.emplace_back(input[i].size());
  }
  return output;
}

struct ProblemUtilities_StringLengthsBackwards { 
  using this_t = ProblemUtilities_StringLengthsBackwards;
  using input_t = Problem_StringLengthsBackwards_input_t;
  using output_t = Problem_StringLengthsBackwards_output_t;
  
  using testcase_set_t = TestCaseSet<input_t, output_t>;
  
  testcase_set_t testing_set;
  testcase_set_t training_set;

  // --- Useful during a test evaluation ---
  bool submitted;
  emp::vector<size_t> submitted_val;
  
  ProblemUtilities_StringLengthsBackwards()
    : testing_set(this_t::LoadTestCaseFromLine),
      training_set(this_t::LoadTestCaseFromLine),
      submitted(false), submitted_val()
  { ; }

  ~ProblemUtilities_StringLengthsBackwards() { ; }

  testcase_set_t & GetTestingSet() { return testing_set; }
  testcase_set_t & GetTrainingSet() { return training_set; }

  input_t & GetTestInput(bool training, size_t id) { 
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetInput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetInput(id); }
  }
  output_t & GetTestOutput(bool training, size_t id) {
    if (training) { emp_assert(id < training_set.GetSize()); return training_set.GetOutput(id); }
    else { emp_assert(id < testing_set.GetSize()); return testing_set.GetOutput(id); }
  }

  void ResetTestEval() {
    submitted = false;
    submitted_val.clear();
  }

  void Submit(size_t val) {
    submitted = true;
    submitted_val.emplace_back(val);
  }

  void Submit(const emp::vector<size_t> & vec) {
    submitted = true;
    submitted_val = vec;
  }
  
  static std::pair<input_t, output_t> LoadTestCaseFromLine(const emp::vector<std::string> & line) {
    input_t input;   
    output_t output; 
    
    // Load input.
    // std::cout << "==== LOAD TEST CASE FROM LINE ====" << std::endl;
    // Parse line[0]
    std::string input_str = line[0];
    if (input_str.front() == '[') { input_str.erase(0, 1); }
    if (input_str.back() == ']') { input_str.pop_back(); }
    input_str += "\n";
    
    // std::cout << "Input str =" << input_str << std::endl;
    // std::cout << "Line[0] = " << line[0] << std::endl;
    // std::cout << "(1) Replace commas" << std::endl;
    input_str = StrReplace(input_str, ",", "{COMMA}");                  // Get rid of commas to avoid confusing CSV parser
    // std::cout << "(1.25) Replace escaped backslashes" << std::endl;
    input_str = StrReplace(input_str, "\\\\", "{BSLASH}");              // Get rid of backslashes so we can get rid of escaped quotes to avoid confusing the CSV parser
    // std::cout << "(1.5) Replace escaped uotes" << std::endl; 
    input_str = StrReplace(input_str, "\\\"", "{QUOTE}");               // Get rid of quotes to avoid confusing the parser
    // std::cout << "  Input str = " << input_str << std::endl; 
    
    // We use a csv parser to tackle the input line (after a bit of processing things that confuse the parser...)
    std::istringstream instr(input_str);
    aria::csv::CsvParser parser(instr);
    parser.delimiter(' ');
    parser.quote('"');
    parser.terminator('\n');

    size_t cnt = 0;
    for (auto & row : parser) {
      ++cnt;
      for (auto & field : row) {
        std::string fstr = StrReplace(field, "{COMMA}", ",");
        fstr = StrReplace(fstr, "\\t", "\t");
        fstr = StrReplace(fstr, "\\n", "\n");
        fstr = StrReplace(fstr, "{QUOTE}", "\"");
        fstr = StrReplace(fstr, "{BSLASH}", "\\");
        // std::cout << "- Instr Field=" << fstr << std::endl;
        input.emplace_back(fstr);
      }
    }

    emp::vector<size_t> correct_out;
    emp::vector<size_t> read_out;

    // What do we *expect* the correct output to be?
    for (int i = input.size()-1; i >= 0; --i) { correct_out.emplace_back(input[i].size()); }
    
    // Handle output
    if (line.size() > 1) {
      // What do we read the correct output as?
      emp::vector<std::string> sliced_line = emp::slice(line[1], '\n');
      for (size_t i = 0; i < sliced_line.size(); ++i) {
        read_out.emplace_back(std::atoi(sliced_line[i].c_str()));
      }
      // If read out and correct out are not the same size, throw an error.
      if (read_out.size() != correct_out.size()) {
        std::cout << "ERROR! Read len(" << read_out.size() << ") != gen len (" << correct_out.size() << "). Exiting" << std::endl;
        exit(-1);
      }
      // If read out and correct out are not identical, throw an error.
      if (read_out != correct_out) {
        std::cout << "ERROR! Read output different from calculated output!" << std::endl;
        
        std::cout << "Read out: ";
        for (size_t i = 0; i < read_out.size(); ++i) {
          if (i) std::cout << ",";
          std::cout << read_out[i];
        } std::cout << std::endl;

        std::cout << "Calculated out: ";
        for (size_t i = 0; i < correct_out.size(); ++i) {
          if (i) std::cout << ",";
          std::cout << correct_out[i];
        } std::cout << std::endl;
        
        std::cout << "Exiting." << std::endl;
        exit(-1);
      }

    }

    output = correct_out;

    // std::cout << "output = [";
    // for (size_t i = 0; i < output.size(); ++i) {
      // if (i) std::cout << ",";
      // std::cout << output[i];
    // } std::cout << "]" << std::endl;

    return {input, output};
  }

  std::pair<double, bool> CalcScorePassFail(const output_t & correct_test_output, const output_t & sub) {
    const bool pass = (sub == correct_test_output);
    return {(double)pass, pass};
  }

  std::pair<double, bool> CalcScoreGradient(const output_t & correct_test_output, const output_t & sub) {
    const double max_dist = emp::Max(correct_test_output.size(), sub.size());
    double dist = emp::calc_edit_distance(correct_test_output, sub);
    if (dist == 0) {
      return {1.0, true};
    } else {
      return {(max_dist - dist)/max_dist, false};
    }
  } 

  void PrintTestCSV(std::ostream & os, const input_t & in) const {
    os << "\"";
    os << "[";
    for (size_t i = 0; i < in.size(); ++i) {
      if (i) os << ",";
      os << "{BEGIN-STR}" << in[i] << "{END-STR}";
    }
    os << "]"; 
    os << "\"";
  }
};

#endif