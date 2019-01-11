# Hand-coded solutions for programming synthesis benchmark problems

<!-- TOC -->

- [Problem - Number IO](#problem---number-io)
  - [Number IO - Tag-based arguments](#number-io---tag-based-arguments)
  - [Number IO - Numeric arguments](#number-io---numeric-arguments)
- [Small or Large](#small-or-large)
  - [Small or Large - Tag-based arguments](#small-or-large---tag-based-arguments)
  - [Small or Large - Numeric arguments](#small-or-large---numeric-arguments)
- [For Loop Index](#for-loop-index)
  - [For Loop Index - Tag-based arguments](#for-loop-index---tag-based-arguments)
  - [For Loop Index - Numeric arguments](#for-loop-index---numeric-arguments)
- [Grade](#grade)
  - [Grade - Tag-based arguments](#grade---tag-based-arguments)
  - [Grade - Numeric arguments](#grade---numeric-arguments)
- [Median](#median)
  - [Median - Tag-based arguments](#median---tag-based-arguments)
  - [Median - Numeric arguments](#median---numeric-arguments)
- [Smallest](#smallest)
  - [Smallest - Tag-based arguments](#smallest---tag-based-arguments)
  - [Smallest - Numeric arguments](#smallest---numeric-arguments)

<!-- /TOC -->

## Problem - Number IO

### Number IO - Tag-based arguments

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);
  sol.PushInst("LoadInt-Tag",    {matrix[0], matrix[0], matrix[0]});
  sol.PushInst("LoadDouble-Tag", {matrix[1], matrix[1], matrix[1]});
  sol.PushInst("Add-Tag",        {matrix[0], matrix[1], matrix[2]});
  sol.PushInst("SubmitNum-Tag",  {matrix[2], matrix[2], matrix[2]});
  prog_world->Inject(sol, PROG_POP_SIZE);
```

### Number IO - Numeric arguments

```{C++}
  hardware_t::Program sol(inst_lib);
  sol.PushInst("LoadInt-Num",    {0, 0, 0});
  sol.PushInst("LoadDouble-Num", {1, 1, 1});
  sol.PushInst("Add-Num",        {0, 1, 2});
  sol.PushInst("SubmitNum-Num",  {2, 2, 2});
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Small or Large

### Small or Large - Tag-based arguments

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);
  // Create thresholds to compare n to.
  sol.PushInst("Set-10-Tag",      {matrix[0], matrix[4], matrix[4]});
  sol.PushInst("Set-2-Tag",       {matrix[1], matrix[4], matrix[4]});
  sol.PushInst("Mult-Tag",        {matrix[0], matrix[0], matrix[2]});
  sol.PushInst("Mult-Tag",        {matrix[0], matrix[2], matrix[0]});
  sol.PushInst("Mult-Tag",        {matrix[0], matrix[1], matrix[1]});
  // Load input
  sol.PushInst("LoadInt-Tag",     {matrix[2], matrix[4], matrix[4]});
  // Check if n < 1000
  sol.PushInst("TestNumLess-Tag", {matrix[2], matrix[0], matrix[3]});
  sol.PushInst("If-Tag",          {matrix[3], matrix[4], matrix[4]});
  sol.PushInst("SubmitSmall",     {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("Return",          {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("Close",           {matrix[4], matrix[4], matrix[4]});
  // Check if n < 2000
  sol.PushInst("TestNumLess-Tag", {matrix[2], matrix[1], matrix[3]});
  sol.PushInst("If-Tag",          {matrix[3], matrix[4], matrix[4]});
  sol.PushInst("SubmitNone",      {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("Return",          {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("Close",           {matrix[4], matrix[4], matrix[4]});
  // n must be >= 2000
  sol.PushInst("SubmitLarge",     {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("Return",          {matrix[4], matrix[4], matrix[4]});
  prog_world->Inject(sol, PROG_POP_SIZE);
```

### Small or Large - Numeric arguments

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);
  // Create thresholds to compare n to.
  sol.PushInst("Set-10-Num",      {0, 4, 4});
  sol.PushInst("Set-2-Num",       {1, 4, 4});
  sol.PushInst("Mult-Num",        {0, 0, 2});
  sol.PushInst("Mult-Num",        {0, 2, 0});
  sol.PushInst("Mult-Num",        {0, 1, 1});
  // Load input
  sol.PushInst("LoadInt-Num",     {2, 4, 4});
  // Check if n < 1000
  sol.PushInst("TestNumLess-Num", {2, 0, 3});
  sol.PushInst("If-Num",          {3, 4, 4});
  sol.PushInst("SubmitSmall",     {4, 4, 4});
  sol.PushInst("Return",          {4, 4, 4});
  sol.PushInst("Close",           {4, 4, 4});
  // Check if n < 2000
  sol.PushInst("TestNumLess-Num", {2, 1, 3});
  sol.PushInst("If-Num",          {3, 4, 4});
  sol.PushInst("SubmitNone",      {4, 4, 4});
  sol.PushInst("Return",          {4, 4, 4});
  sol.PushInst("Close",           {4, 4, 4});
  // n must be >= 2000
  sol.PushInst("SubmitLarge",     {4, 4, 4});
  sol.PushInst("Return",          {4, 4, 4});
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## For Loop Index

### For Loop Index - Tag-based arguments

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("CopyMem-Tag",     {matrix[0], matrix[4], matrix[7]});
  sol.PushInst("Inc-Tag",         {matrix[5], matrix[7], matrix[7]});
  sol.PushInst("While-Tag",       {matrix[5], matrix[7], matrix[7]});
  sol.PushInst("SubmitNum-Tag",   {matrix[4], matrix[7], matrix[7]});
  sol.PushInst("Add-Tag",         {matrix[4], matrix[2], matrix[4]});
  sol.PushInst("TestNumLess-Tag", {matrix[4], matrix[1], matrix[5]});
  sol.PushInst("Close",           {matrix[7], matrix[7], matrix[7]});

  prog_world->Inject(sol, PROG_POP_SIZE);
```

### For Loop Index - Numeric arguments

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("CopyMem-Num",     {0, 4, 7});
  sol.PushInst("Inc-Num",         {5, 7, 7});
  sol.PushInst("While-Num",       {5, 7, 7});
  sol.PushInst("SubmitNum-Num",   {4, 7, 7});
  sol.PushInst("Add-Num",         {4, 2, 4});
  sol.PushInst("TestNumLess-Num", {4, 1, 5});
  sol.PushInst("Close",           {7, 7, 7});

  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Grade

### Grade - Tag-based arguments

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);
  sol.PushInst("LoadThreshA-Tag",        {matrix[0], matrix[8], matrix[8]});
  sol.PushInst("LoadThreshB-Tag",        {matrix[1], matrix[8], matrix[8]});
  sol.PushInst("LoadThreshC-Tag",        {matrix[2], matrix[8], matrix[8]});
  sol.PushInst("LoadThreshD-Tag",        {matrix[3], matrix[8], matrix[8]});
  sol.PushInst("LoadGrade-Tag",              {matrix[4], matrix[8], matrix[8]});
  sol.PushInst("TestNumGreaterTEqu-Tag", {matrix[4], matrix[0], matrix[5]});
  sol.PushInst("If-Tag",                 {matrix[5], matrix[8], matrix[8]});
    sol.PushInst("SubmitA",              {matrix[8], matrix[8], matrix[8]});
    sol.PushInst("Return",               {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("Close",                  {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("TestNumGreaterTEqu-Tag", {matrix[4], matrix[1], matrix[5]});
  sol.PushInst("If-Tag",                 {matrix[5], matrix[8], matrix[8]});
    sol.PushInst("SubmitB",              {matrix[8], matrix[8], matrix[8]});
    sol.PushInst("Return",               {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("Close",                  {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("TestNumGreaterTEqu-Tag", {matrix[4], matrix[2], matrix[5]});
  sol.PushInst("If-Tag",                 {matrix[5], matrix[8], matrix[8]});
    sol.PushInst("SubmitC",              {matrix[8], matrix[8], matrix[8]});
    sol.PushInst("Return",               {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("Close",                  {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("TestNumGreaterTEqu-Tag", {matrix[4], matrix[3], matrix[5]});
  sol.PushInst("If-Tag",                 {matrix[5], matrix[8], matrix[8]});
    sol.PushInst("SubmitD",              {matrix[8], matrix[8], matrix[8]});
    sol.PushInst("Return",               {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("Close",                  {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("SubmitF",                {matrix[8], matrix[8], matrix[8]});
  
  prog_world->Inject(sol, PROG_POP_SIZE);
```

### Grade - Numeric arguments

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);
  sol.PushInst("LoadThreshA-Num",        {0, 8, 8});
  sol.PushInst("LoadThreshB-Num",        {1, 8, 8});
  sol.PushInst("LoadThreshC-Num",        {2, 8, 8});
  sol.PushInst("LoadThreshD-Num",        {3, 8, 8});
  sol.PushInst("LoadGrade-Num",              {4, 8, 8});
  sol.PushInst("TestNumGreaterTEqu-Num", {4, 0, 5});
  sol.PushInst("If-Num",                 {5, 8, 8});
    sol.PushInst("SubmitA",              {8, 8, 8});
    sol.PushInst("Return",               {8, 8, 8});
  sol.PushInst("Close",                  {8, 8, 8});
  sol.PushInst("TestNumGreaterTEqu-Num", {4, 1, 5});
  sol.PushInst("If-Num",                 {5, 8, 8});
    sol.PushInst("SubmitB",              {8, 8, 8});
    sol.PushInst("Return",               {8, 8, 8});
  sol.PushInst("Close",                  {8, 8, 8});
  sol.PushInst("TestNumGreaterTEqu-Num", {4, 2, 5});
  sol.PushInst("If-Num",                 {5, 8, 8});
    sol.PushInst("SubmitC",              {8, 8, 8});
    sol.PushInst("Return",               {8, 8, 8});
  sol.PushInst("Close",                  {8, 8, 8});
  sol.PushInst("TestNumGreaterTEqu-Num", {4, 3, 5});
  sol.PushInst("If-Num",                 {5, 8, 8});
    sol.PushInst("SubmitD",              {8, 8, 8});
    sol.PushInst("Return",               {8, 8, 8});
  sol.PushInst("Close",                  {8, 8, 8});
  sol.PushInst("SubmitF",                {8, 8, 8});
  
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Median

### Median - Tag-based arguments

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  // todo

  prog_world->Inject(sol, PROG_POP_SIZE);
```

### Median - Numeric arguments

## Smallest

### Smallest - Tag-based arguments

### Smallest - Numeric arguments

