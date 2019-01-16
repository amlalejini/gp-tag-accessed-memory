#ifndef PROGRAM_SYNTHESIS_CONFIG_H
#define PROGRAM_SYNTHESIS_CONFIG_H

#include "config/config.h"

EMP_BUILD_CONFIG(ProgramSynthesisConfig, 
  GROUP(DEFAULT_GROUP, "General settings"),
  VALUE(SEED, int, 0, "Random number seed (-1 for based on time)"),
  VALUE(GENERATIONS, size_t, 10000, "How many generations should we run for?"),
  VALUE(PROGRAM_ARGUMENT_MODE, size_t, 0, "What type of arguments should programs make use of? \n0: Tag arguments \n1: Numeric arguments \n2: Both tag and numeric arguments"),
  VALUE(PROGRAM_ARGUMENTS_TYPE_SEARCH, bool, false, "Should instruction arguments search memory for the closest matching position OF THE CORRECT TYPE?"),
  VALUE(PROG_POP_SIZE, size_t, 1024, "Population size for programs"),
  VALUE(PROBLEM, std::string, "number-io", "Which problem to use?"),
  VALUE(BENCHMARK_DATA_DIR, std::string, "../data/prog-synth-examples", "Location to look for problem test case data."),

  GROUP(PROGRAM_GROUP, "General settings specific to programs."),
  VALUE(USE_MODULES, bool, false, "Allow modules?"),
  VALUE(MIN_PROG_SIZE, size_t, 1, "Minimum program size"),
  VALUE(MAX_PROG_SIZE, size_t, 128, "Maximum program size"),
  VALUE(PROG_EVAL_TIME, size_t, 256, "How many clock cycles should we give a program during a test?"),
  VALUE(PROG_MUT__PER_BIT_FLIP, double, 0.001, "Program per-bit flip rate."),
  VALUE(PROG_MUT__PER_NUMERIC_ARG_SUB, double, 0.001, "Program numeric argument substitution rate."),
  VALUE(PROG_MUT__PER_INST_SUB, double, 0.005, "Program per-instruction substitution mutation rate."),
  VALUE(PROG_MUT__PER_INST_INS, double, 0.005, "Program per-instruction insertion mutation rate."),
  VALUE(PROG_MUT__PER_INST_DEL, double, 0.005, "Program per-instruction deletion mutation rate."),
  VALUE(PROG_MUT__PER_PROG_SLIP, double, 0.05, "Program per-program slip mutation rate."),
  VALUE(PROG_MUT__PER_MOD_DUP, double, 0.05, "Program per-module whole-module duplication rate."),
  VALUE(PROG_MUT__PER_MOD_DEL, double, 0.05, "Program per-module whole-module deletion rate."),

  GROUP(HARDWARE_GROUP, "Settings specific to TagLGP virtual hardware"),
  VALUE(MIN_TAG_SPECIFICITY, double, 0.0, "What is the minimum tag similarity required for a tag to successfully reference another tag?"),
  VALUE(MAX_CALL_DEPTH, size_t, 128, "Maximum depth of hardware's call stack."),

  GROUP(SELECTION_GROUP, "Setting specific to selection"),
  VALUE(LEXICASE_MAX_FUNS, size_t, 0, "Max functions for lexicase selection"),

  GROUP(DATA_COLLECTION_GROUP, "Settings specific to data collection."),
  VALUE(DATA_DIRECTORY, std::string, "./output", "Where should we dump output files?"),
  VALUE(SNAPSHOT_INTERVAL, size_t, 1000, "How often should we take population snapshots?"),
  VALUE(SUMMARY_STATS_INTERVAL, size_t, 1000, "How often should we output summary stats?"),
)

#endif