#ifndef TAG_LINEAR_GP_UTILITIES_H
#define TAG_LINEAR_GP_UTILITIES_H

#include "tools/Random.h"
#include "TagLinearGP.h"
#include "Utilities.h"

namespace TagLGP {

  /// ====== Generators ======

  /// Utility function to generate a random TagGP program with TAG ARGUMENTS.
  template<size_t TAG_WIDTH, size_t MAX_ARG_CNT>
  typename TagLinearGP_TW<TAG_WIDTH>::Program GenRandTagGPProgram(emp::Random & rnd,
                                                                  emp::Ptr<InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT>> inst_lib,
                                                                  size_t MIN_LEN=1, size_t MAX_LEN=100) 
  {
    emp_assert(MIN_LEN < MAX_LEN);
    emp_assert(inst_lib->GetSize());
    using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
    using program_t = typename hardware_t::program_t;
    
    program_t prog(inst_lib);
    size_t N = rnd.GetUInt(MIN_LEN, MAX_LEN+1);
    for (size_t i = 0; i < N; ++i) {
      prog.PushInst(GenRandTagGPInst(rnd, *inst_lib));
    }
    return prog;
  }

  /// Utility function that generates a random TagGP instruction with TAG ARGUMENTS.
  template<size_t TAG_WIDTH, size_t MAX_ARG_CNT> 
  typename TagLinearGP_TW<TAG_WIDTH>::Instruction GenRandTagGPInst(emp::Random & rnd, 
                                                                   const InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT> & inst_lib) 
  {
    emp_assert(inst_lib.GetSize() > 0, "Instruction library must have at least one instruction definition before being used to generate a random instruction.");
    
    // using inst_lib_t = InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT>;
    using hardware_t = TagLinearGP_TW<TAG_WIDTH>;
    using inst_t = typename hardware_t::inst_t;
    using tag_t = typename hardware_t::tag_t;

    // Generate tag arguments.
    emp::vector<tag_t> tags;
    for (size_t i = 0; i < MAX_ARG_CNT; ++i) {
      tags.emplace_back(rnd, 0.5);
    }
    
    return inst_t(rnd.GetUInt(inst_lib.GetSize()), tags);
  }

  /// Utility function that generates a random TagGP instruction with NUMERIC ARGUMENTS.
  template<size_t TAG_WIDTH, size_t MAX_ARG_CNT> 
  typename TagLinearGP_TW<TAG_WIDTH>::Instruction GenRandTagGPInst_NumArgs(emp::Random & rnd, 
                                                                           const InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT> & inst_lib,
                                                                           size_t max_arg) 
  {
    emp_assert(inst_lib.GetSize() > 0, "Instruction library must have at least one instruction definition before being used to generate a random instruction.");
    
    // using inst_lib_t = InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT>;
    using hardware_t = TagLinearGP_TW<TAG_WIDTH>;
    using inst_t = typename hardware_t::inst_t;

    // Generate numeric arguments.
    emp::vector<size_t> args;
    for (size_t i = 0; i < MAX_ARG_CNT; ++i) {
      args.emplace_back(rnd.GetUInt(0, max_arg+1));
    }

    return inst_t(rnd.GetUInt(inst_lib.GetSize()), args);
  }

  /// Utility function that generates a random TagGP program with NUMERIC ARGUMENTS.
  template<size_t TAG_WIDTH, size_t MAX_ARG_CNT>
  typename TagLinearGP_TW<TAG_WIDTH>::Program GenRandTagGPProgram_NumArgs(emp::Random & rnd,
                                                                          emp::Ptr<InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT>> inst_lib,
                                                                          size_t max_arg,
                                                                          size_t MIN_LEN=1, size_t MAX_LEN=100) 
  {
    emp_assert(MIN_LEN < MAX_LEN);
    emp_assert(inst_lib->GetSize());
    using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
    using program_t = typename hardware_t::program_t;
    
    program_t prog(inst_lib);
    size_t N = rnd.GetUInt(MIN_LEN, MAX_LEN+1);
    for (size_t i = 0; i < N; ++i) {
      prog.PushInst(GenRandTagGPInst_NumArgs(rnd, *inst_lib, max_arg));
    }
    return prog;
  }

  /// Utility function that generates a random TagGP instruction with NUMERIC AND TAG ARGUMENTS.
  template<size_t TAG_WIDTH, size_t MAX_ARG_CNT> 
  typename TagLinearGP_TW<TAG_WIDTH>::Instruction GenRandTagGPInst_TagAndNumArgs(emp::Random & rnd, 
                                                                                 const InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT> & inst_lib,
                                                                                 size_t max_arg) 
  {
    emp_assert(inst_lib.GetSize() > 0, "Instruction library must have at least one instruction definition before being used to generate a random instruction.");
    
    // using inst_lib_t = InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT>;
    using hardware_t = TagLinearGP_TW<TAG_WIDTH>;
    using inst_t = typename hardware_t::inst_t;
    using tag_t = typename hardware_t::tag_t;

    // Generate numeric arguments.
    emp::vector<size_t> args;
    for (size_t i = 0; i < MAX_ARG_CNT; ++i) {
      args.emplace_back(rnd.GetUInt(0, max_arg+1));
    }

    // Generate tag arguments.
    emp::vector<tag_t> tags;
    for (size_t i = 0; i < MAX_ARG_CNT; ++i) {
      tags.emplace_back(rnd, 0.5);
    }

    return inst_t(rnd.GetUInt(inst_lib.GetSize()), tags, args);
  }

  /// Utility function that generates a random TagGP program with NUMERIC AND TAG ARGUMENTS.
  template<size_t TAG_WIDTH, size_t MAX_ARG_CNT>
  typename TagLinearGP_TW<TAG_WIDTH>::Program GenRandTagGPProgram_TagAndNumArgs(emp::Random & rnd,
                                                                                emp::Ptr<InstLib<TagLinearGP_TW<TAG_WIDTH>, MAX_ARG_CNT>> inst_lib,
                                                                                size_t max_arg,
                                                                                size_t MIN_LEN=1, size_t MAX_LEN=100) 
  {
    emp_assert(MIN_LEN < MAX_LEN);
    emp_assert(inst_lib->GetSize());
    using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>;
    using program_t = typename hardware_t::program_t;
    
    program_t prog(inst_lib);
    size_t N = rnd.GetUInt(MIN_LEN, MAX_LEN+1);
    for (size_t i = 0; i < N; ++i) {
      prog.PushInst(GenRandTagGPInst_TagAndNumArgs(rnd, *inst_lib, max_arg));
    }
    return prog;
  }

}

#endif