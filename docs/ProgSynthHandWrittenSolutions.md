# Hand-coded solutions for programming synthesis benchmark problems

<!-- TOC -->

- [Problem - Number IO](#problem---number-io)
  - [Tag-based arguments](#tag-based-arguments)
  - [Numeric arguments](#numeric-arguments)

<!-- /TOC -->

## Problem - Number IO

### Tag-based arguments

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);
  sol.PushInst("LoadInt-Tag",    {matrix[0], matrix[0], matrix[0]});
  sol.PushInst("LoadDouble-Tag", {matrix[1], matrix[1], matrix[1]});
  sol.PushInst("Add-Tag",        {matrix[0], matrix[1], matrix[2]});
  sol.PushInst("SubmitNum-Tag",  {matrix[2], matrix[2], matrix[2]});
  prog_world->Inject(sol, PROG_POP_SIZE);
```

### Numeric arguments

```{C++}
  hardware_t::Program sol(inst_lib);
  sol.PushInst("LoadInt-Num",    {0, 0, 0});
  sol.PushInst("LoadDouble-Num", {1, 1, 1});
  sol.PushInst("Add-Num",        {0, 1, 2});
  sol.PushInst("SubmitNum-Num",  {2, 2, 2});
  prog_world->Inject(sol, PROG_POP_SIZE);
```
