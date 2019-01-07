#ifndef TAGMEM_MUTATORS_H
#define TAGMEM_MUTATORS_H

#include "base/vector.h"
#include "tools/Random.h"
#include "tools/random_utils.h"

#include "TagLinearGP.h"
#include "TagLinearGP_Utilities.h"

/// Mutator for tag LGP programs
template<size_t TAG_WIDTH>
struct TagLGPMutator {
  using hardware_t = TagLGP::TagLinearGP_TW<TAG_WIDTH>; 
  using program_t = typename hardware_t::Program;
  using tag_t = typename hardware_t::tag_t;
  using inst_t = typename hardware_t::inst_t;
  using inst_lib_t = typename hardware_t::inst_lib_t;
  using inst_prop_t = typename hardware_t::inst_prop_t;
  
  size_t MAX_PROGRAM_LEN;
  size_t MIN_PROGRAM_LEN;
  size_t MAX_NUMERIC_ARG;

  bool USE_NUMERIC_ARGUMENTS;
  bool USE_TAG_ARGUMENTS;

  double PER_BIT_FLIP;
  double PER_NUMERIC_ARG_SUB;
  double PER_INST_SUB;
  double PER_INST_INS;
  double PER_INST_DEL;
  double PER_PROG_SLIP;
  double PER_MOD_DUP;
  double PER_MOD_DEL;
  
  enum ModuleMutType { DUP=0, DEL, NONE };

  struct ModuleInfo {
    int begin;  // Includes ModDef instruction
    int end;
    ModuleMutType mut;
    emp::vector<size_t> positions;
    ModuleInfo(int _b, int _e) : begin(_b), end(_e), mut(ModuleMutType::NONE), positions() { ; }
    size_t GetSize() const { return positions.size(); }
  };

  size_t Mutate(emp::Random & rnd, program_t & program) {
    const inst_lib_t & ilib = program.GetInstLib();
    
    int expected_size = program.GetSize();
    size_t mut_cnt = 0;
    int rhead = 0;

    // assert that expected size must be in correct range
    /////////////////////////////////////////// 
    // Whole module duplications/deletions
    if (PER_MOD_DEL != 0.0 || PER_MOD_DUP != 0.0) {
      
      // Scan for modules
      size_t mod_mut_cnt = 0;
      emp::vector<ModuleInfo> modules;  // Modules we're going to dup or delete.
      emp::vector<size_t> module_membership(program.GetSize(), 0);
      emp::vector<bool> dels(program.GetSize(), false);
      emp::vector<size_t> danglers;
      
      for (size_t i = 0; i < program.GetSize(); ++i) {
        // inst_t & inst = program[i];
        if (ilib.HasProperty(program[i].id, inst_prop_t::MODULE)) {
          if (modules.size()) { modules.back().end = i; }
          const size_t mod_id = modules.size();
          modules.emplace_back(i, -1);
          modules.back().positions.emplace_back(i); // Add start position to module
          module_membership[i] = mod_id;
        } else {
          // Not a new module definition
          if (modules.size()) { modules.back().positions.emplace_back(i); module_membership[i] = modules.size()-1; }
          else { danglers.emplace_back(i); }
        }
      }

      // Take care of danglers.
      if (modules.size()) {
        if (modules[0].begin == 0) { modules.back().end = program.GetSize(); }  // Case where last module does not wrap.
        else { modules.back().end = modules[0].begin-1; }                       // Last module wraps around.
        for (size_t i = 0; i < danglers.size(); ++i) {
          modules.back().positions.emplace_back(danglers[i]); 
          module_membership[i] = modules.size()-1;
        }
      }

      // Did we do the above thing correctly?
      // std::cout << "Modules found: " << modules.size() << std::endl;
      // for (size_t i = 0; i < modules.size(); ++i) {
      //   std::cout << "  Module " << i << "["<<modules[i].begin << "," << modules[i].end << ")" << std::endl;
      //   for (size_t p = 0; p < modules[i].positions.size(); ++p) {
      //     if (p == 0) std::cout << "    => [";
      //     else std::cout << ",";
      //     std::cout << modules[i].positions[p];
      //   }
      //   std::cout << "]" << std::endl;
      // }

      // Are we mutating any of the modules?
      for (size_t i = 0; i < modules.size(); ++i) {
        bool mod_dup = rnd.P(PER_MOD_DUP);
        bool mod_del = rnd.P(PER_MOD_DEL);
        if (mod_dup && mod_del) { mod_dup = false; mod_dup = false; } // If we would both dup and delete module, do nothing instead.
        if (mod_dup && ((expected_size + modules[i].GetSize()) <= MAX_PROGRAM_LEN) ) { 
          mod_mut_cnt++; 
          mut_cnt++; 
          modules[i].mut = ModuleMutType::DUP; 
          expected_size += modules[i].GetSize(); 
        }
        if (mod_del && (((int)expected_size - (int)modules[i].GetSize()) >= (int)MIN_PROGRAM_LEN) ) { 
          mod_mut_cnt++; 
          mut_cnt++; 
          modules[i].mut = ModuleMutType::DEL; 
          expected_size -= modules[i].GetSize(); 
          for (size_t p = 0; p < modules[i].GetSize(); ++p) dels[modules[i].positions[p]] = true;
        }
      }
      
      // Do all of the dups/deletions
      if (mod_mut_cnt > 0) {
        program_t new_program(program.GetInstLibPtr()); // Program we're writing to. (will be copied over.) 
        // for (rhead = 0; rhead < program.GetSize(); ++rhead) {
        rhead = 0;
        while ((int)new_program.GetSize() < expected_size) {
          const size_t cur_module = module_membership[rhead];
          // Did we find a new module?
          if (ilib.HasProperty(program[rhead].id, inst_prop_t::MODULE)) {
            // Should we duplicate cur_module?
            if (modules[cur_module].mut == ModuleMutType::DUP) {
              // We're supposed to duplicate current module.
              for (size_t i = 0; i < modules[cur_module].GetSize(); ++i) {
                new_program.PushInst(program[modules[cur_module].positions[i]]);
              }
            }
          }
          if (!dels[rhead]) {
            new_program.PushInst(program[rhead]);
          }
          ++rhead;
        }
        program = new_program;
      }

    }

    // Slip? -> If so, where?
    bool slip = false;
    bool slip_dup = false;
    bool slip_del = false;
    size_t slip_begin=0;
    size_t slip_end=0;
    // int slip_dup_size = 0;
    // int slip_del_size = 0;
    if (rnd.P(PER_PROG_SLIP)) {
      slip_begin = rnd.GetUInt(program.GetSize());
      slip_end = rnd.GetUInt(program.GetSize());
      slip_dup = slip_begin < slip_end;
      slip_del = slip_begin > slip_end;
      slip = slip_dup || slip_del; // If we may dup or del, we're slipping! (well, maybe - need to check that slip would't break constraints)
      if (slip) {
        const int slip_dup_size = 1 + (int)slip_end - (int)slip_begin;
        const int slip_del_size = 1 + (int)slip_begin - (int)slip_end;

        if (slip_dup && ((expected_size+slip_dup_size) > (int)MAX_PROGRAM_LEN)) { 
          // If slip-dup would break constraints, don't slip.
          slip = false; slip_dup=false; 
        } 
        if (slip_dup) { expected_size += slip_dup_size; }

        if (slip_del && ((expected_size-slip_del_size) < (int)MIN_PROGRAM_LEN)) { 
          // If slip-del would break constraints, don't slip.
          slip = false; slip_del = false; 
        } 
        if (slip_del) { expected_size -= slip_del_size; }

      }
    }

    program_t new_program(program.GetInstLibPtr());
    for (rhead = 0; rhead < (int)program.GetSize(); ++rhead) {
      // Check for slip.
      if (slip_dup && rhead == (int)slip_end) {
        // Copy over inst.
        new_program.PushInst(program[rhead]);
        // Duplicate slip segment.
        for (size_t si = slip_begin; si <= slip_end; ++si) {
          new_program.PushInst(program[si]);
        }
        mut_cnt++;
        continue;
      } else if (slip_del && rhead == (int)slip_end) {
        mut_cnt++;
        rhead = slip_begin;
        continue;
      }
      
      // Instruction deletion
      if (rnd.P(PER_INST_DEL) && ((expected_size-1)>=(int)MIN_PROGRAM_LEN)) {
        --expected_size;
        ++mut_cnt;
        continue;
      }

      // Copy instruction at read head
      const size_t whead = new_program.GetSize();
      new_program.PushInst(program[rhead]);

      // Instruction insertion
      if (rnd.P(PER_INST_INS) && ((expected_size+1)<=(int)MAX_PROGRAM_LEN)) {
        ++expected_size;
        ++mut_cnt;

        if (USE_NUMERIC_ARGUMENTS && !USE_TAG_ARGUMENTS) {
          new_program.PushInst(TagLGP::GenRandTagGPInst_NumArgs(rnd, ilib, MAX_NUMERIC_ARG));
        } else if (!USE_NUMERIC_ARGUMENTS && USE_TAG_ARGUMENTS) {
          new_program.PushInst(TagLGP::GenRandTagGPInst(rnd, ilib));
        } else { // Allowing for BOTH types of arguments
          new_program.PushInst(TagLGP::GenRandTagGPInst_TagAndNumArgs(rnd, ilib, MAX_NUMERIC_ARG));
        }
      }

      // Instruction substitution
      if (rnd.P(PER_INST_SUB)) {
        ++mut_cnt;
        new_program[whead].id = rnd.GetUInt(ilib.GetSize());
      }

      // Argument substitutions
      // - Numeric arguments (if there are none, will never enter for loop)
      for (size_t arg = 0; arg < new_program[whead].arg_nums.size(); ++arg) {
        if (rnd.P(PER_NUMERIC_ARG_SUB)) {
          new_program[whead].arg_nums[arg] = rnd.GetUInt(0, MAX_NUMERIC_ARG+1);
        }
      }
      // - tag arguments (if there are none, will never enter for loop)
      for (size_t arg = 0; arg < new_program[whead].arg_tags.size(); ++arg) {
        tag_t & tag = new_program[whead].arg_tags[arg];
        for (size_t k = 0; k < tag.GetSize(); ++k) {
          if (rnd.P(PER_BIT_FLIP)) {
            ++mut_cnt;
            tag.Toggle(k);
          }
        }
      }
    }

    program = new_program;
    return mut_cnt;
  }



};

#endif