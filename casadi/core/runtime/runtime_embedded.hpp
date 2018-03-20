namespace casadi {
const char * codegen_str_copy =
"void CASADI_PREFIX(copy)(const real_t* x, int n, real_t* y) {\n"
"  int i;\n"
"  if (y) {\n"
"    if (x) {\n"
"      for (i=0; i<n; ++i) *y++ = *x++;\n"
"    } else {\n"
"      for (i=0; i<n; ++i) *y++ = 0.;\n"
"    }\n"
"  }\n"
"}\n";

const char * codegen_str_swap =
"void CASADI_PREFIX(swap)(int n, real_t* x, int inc_x, real_t* y, int inc_y) {\n"
"  real_t t;\n"
"  int i;\n"
"  for (i=0; i<n; ++i) {\n"
"    t = *x;\n"
"    *x = *y;\n"
"    *y = t;\n"
"    x += inc_x;\n"
"    y += inc_y;\n"
"  }\n"
"}\n";

const char * codegen_str_project =
"void CASADI_PREFIX(project)(const real_t* x, const int* sp_x, real_t* y, const int* sp_y, real_t* w) {\n"
"  int ncol_x = sp_x[1];\n"
"  const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;\n"
"  int ncol_y = sp_y[1];\n"
"  const int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;\n"
"  /* Loop over columns of x and y */\n"
"  int i, el;\n"
"  for (i=0; i<ncol_x; ++i) {\n"
"    /* Zero out requested entries in y */\n"
"    for (el=colind_y[i]; el<colind_y[i+1]; ++el) w[row_y[el]] = 0;\n"
"    /* Set x entries */\n"
"    for (el=colind_x[i]; el<colind_x[i+1]; ++el) w[row_x[el]] = x[el];\n"
"    /* Retrieve requested entries in y */\n"
"    for (el=colind_y[i]; el<colind_y[i+1]; ++el) y[el] = w[row_y[el]];\n"
"  }\n"
"}\n";

const char * codegen_str_mproject =
"void CASADI_PREFIX(mproject)(real_t factor, const real_t* x, const int* sp_x, real_t* y, const int* sp_y, real_t* w) {\n"
"  int ncol_x = sp_x[1];\n"
"  const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;\n"
"  int ncol_y = sp_y[1];\n"
"  const int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;\n"
"  /* Loop over columns of x and y */\n"
"  int i, el;\n"
"  for (i=0; i<ncol_x; ++i) {\n"
"    /* Zero out requested entries in y */\n"
"    for (el=colind_y[i]; el<colind_y[i+1]; ++el) w[row_y[el]] = 0;\n"
"    /* Set x entries */\n"
"    for (el=colind_x[i]; el<colind_x[i+1]; ++el) w[row_x[el]] = x[el];\n"
"    /* Retrieve requested entries in y */\n"
"    for (el=colind_y[i]; el<colind_y[i+1]; ++el) y[el] = factor*w[row_y[el]];\n"
"  }\n"
"}\n";

const char * codegen_str_maddproject =
"void CASADI_PREFIX(maddproject)(real_t factor, const real_t* x, const int* sp_x, real_t* y, const int* sp_y, real_t* w) {\n"
"  int ncol_x = sp_x[1];\n"
"  const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;\n"
"  int ncol_y = sp_y[1];\n"
"  const int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;\n"
"  /* Loop over columns of x and y */\n"
"  int i, el;\n"
"  for (i=0; i<ncol_x; ++i) {\n"
"    /* Zero out requested entries in y */\n"
"    for (el=colind_y[i]; el<colind_y[i+1]; ++el) w[row_y[el]] = 0;\n"
"    /* Set x entries */\n"
"    for (el=colind_x[i]; el<colind_x[i+1]; ++el) w[row_x[el]] = x[el];\n"
"    /* Retrieve requested entries in y */\n"
"    for (el=colind_y[i]; el<colind_y[i+1]; ++el) y[el] += factor*w[row_y[el]];\n"
"  }\n"
"}\n";

const char * codegen_str_dense_transfer =
"void CASADI_PREFIX(dense_transfer)(real_t factor, const real_t* x, const int* sp_x, real_t* y, const int* sp_y, real_t* w) {\n"
"  CASADI_PREFIX(sparsify)(x, w, sp_x, false);\n"
"  int nrow_y = sp_y[0];\n"
"  int ncol_y = sp_y[1];\n"
"  const int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;\n"
"  /* Loop over columns of y */\n"
"  int i, el;\n"
"  for (i=0; i<ncol_y; ++i) {\n"
"    for (el=colind_y[i]; el<colind_y[i+1]; ++el) y[nrow_y*i + row_y[el]] += factor*(*w++);\n"
"  }\n"
"}\n";

const char * codegen_str_densify =
"void CASADI_PREFIX(densify)(const real1_t* x, const int* sp_x, real2_t* y, int tr) {\n"
"  /* Quick return - output ignored */\n"
"  if (!y) return;\n"
"  int nrow_x = sp_x[0], ncol_x = sp_x[1];\n"
"  const int *colind_x = sp_x+2, *row_x = sp_x+ncol_x+3;\n"
"  /* Zero out return value */\n"
"  CASADI_PREFIX(fill)(y, nrow_x*ncol_x, CASADI_CAST(real2_t, 0));\n"
"  /* Quick return - input is zero */\n"
"  if (!x) return;\n"
"  /* Copy nonzeros */\n"
"  int i, el;\n"
"  if (tr) {\n"
"    for (i=0; i<ncol_x; ++i) {\n"
"      for (el=colind_x[i]; el!=colind_x[i+1]; ++el) {\n"
"        y[i + row_x[el]*ncol_x] = CASADI_CAST(real2_t, *x++);\n"
"      }\n"
"    }\n"
"  } else {\n"
"    for (i=0; i<ncol_x; ++i) {\n"
"      for (el=colind_x[i]; el!=colind_x[i+1]; ++el) {\n"
"        y[row_x[el] + i*nrow_x] = CASADI_CAST(real2_t, *x++);\n"
"      }\n"
"    }\n"
"  }\n"
"}\n";

const char * codegen_str_sparsify =
"void CASADI_PREFIX(sparsify)(const real1_t* x, real2_t* y, const int* sp_y, int tr) {\n"
"  int nrow_y = sp_y[0], ncol_y = sp_y[1];\n"
"  const int *colind_y = sp_y+2, *row_y = sp_y+ncol_y+3;\n"
"  int i, el;\n"
"  if (tr) {\n"
"    for (i=0; i<ncol_y; ++i) {\n"
"      for (el=colind_y[i]; el!=colind_y[i+1]; ++el) {\n"
"        *y++ = CASADI_CAST(real2_t, x[i + row_y[el]*ncol_y]);\n"
"      }\n"
"    }\n"
"  } else {\n"
"    for (i=0; i<ncol_y; ++i) {\n"
"      for (el=colind_y[i]; el!=colind_y[i+1]; ++el) {\n"
"        *y++ = CASADI_CAST(real2_t, x[row_y[el] + i*nrow_y]);\n"
"      }\n"
"    }\n"
"  }\n"
"}\n";

const char * codegen_str_scal =
"void CASADI_PREFIX(scal)(int n, real_t alpha, real_t* x) {\n"
"  int i;\n"
"  for (i=0; i<n; ++i) *x++ *= alpha;\n"
"}\n";

const char * codegen_str_axpy =
"void CASADI_PREFIX(axpy)(int n, real_t alpha, const real_t* x, real_t* y) {\n"
"  int i;\n"
"  for (i=0; i<n; ++i) *y++ += alpha**x++;\n"
"}\n";

const char * codegen_str_dot =
"real_t CASADI_PREFIX(dot)(int n, const real_t* x, const real_t* y) {\n"
"  real_t r = 0;\n"
"  int i;\n"
"  for (i=0; i<n; ++i) r += *x++ * *y++;\n"
"  return r;\n"
"}\n";

const char * codegen_str_max_viol =
"real_t CASADI_PREFIX(max_viol)(int n, const real_t* x, const real_t* lb, const real_t* ub) {\n"
"  real_t r = 0;\n"
"  const real_t zero = 0;\n"
"  int i;\n"
"  for (i=0; i<n; ++i) {\n"
"    real_t x_i = x ? *x++ : zero;\n"
"    real_t lb_i = lb ? *lb++ : zero;\n"
"    real_t ub_i = ub ? *ub++ : zero;\n"
"    r = fmax(r, fmax(x_i-ub_i, zero));\n"
"    r = fmax(r, fmax(lb_i-x_i, zero));\n"
"  }\n"
"  return r;\n"
"}\n";

const char * codegen_str_sum_viol =
"real_t CASADI_PREFIX(sum_viol)(int n, const real_t* x, const real_t* lb, const real_t* ub) {\n"
"  real_t r = 0;\n"
"  const real_t zero = 0;\n"
"  int i;\n"
"  for (i=0; i<n; ++i) {\n"
"    real_t x_i = x ? *x++ : zero;\n"
"    real_t lb_i = lb ? *lb++ : zero;\n"
"    real_t ub_i = ub ? *ub++ : zero;\n"
"    r += fmax(x_i-ub_i, zero);\n"
"    r += fmax(lb_i-x_i, zero);\n"
"  }\n"
"  return r;\n"
"}\n";

const char * codegen_str_iamax =
"int CASADI_PREFIX(iamax)(int n, const real_t* x, int inc_x) {\n"
"  real_t t;\n"
"  real_t largest_value = -1.0;\n"
"  int largest_index = -1;\n"
"  int i;\n"
"  for (i=0; i<n; ++i) {\n"
"    t = fabs(*x);\n"
"    x += inc_x;\n"
"    if (t>largest_value) {\n"
"      largest_value = t;\n"
"      largest_index = i;\n"
"    }\n"
"  }\n"
"  return largest_index;\n"
"}\n";

const char * codegen_str_fill =
"void CASADI_PREFIX(fill)(real_t* x, int n, real_t alpha) {\n"
"  int i;\n"
"  if (x) {\n"
"    for (i=0; i<n; ++i) *x++ = alpha;\n"
"  }\n"
"}\n";

const char * codegen_str_fill_int =
"void CASADI_PREFIX(fill_int)(int* x, int n, int alpha) {\n"
"  int i;\n"
"  if (x) {\n"
"    for (i=0; i<n; ++i) *x++ = alpha;\n"
"  }\n"
"}\n";

const char * codegen_str_mtimes =
"void CASADI_PREFIX(mtimes)(const real_t* x, const int* sp_x, const real_t* y, const int* sp_y, real_t* z, const int* sp_z, real_t* w, int tr) {\n"
"  /* Get sparsities */\n"
"  int ncol_x = sp_x[1];\n"
"  const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;\n"
"  int ncol_y = sp_y[1];\n"
"  const int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;\n"
"  int ncol_z = sp_z[1];\n"
"  const int *colind_z = sp_z+2, *row_z = sp_z + 2 + ncol_z+1;\n"
"  int cc,kk, kk1;\n"
"  if (tr) {\n"
"    /* Loop over the columns of y and z */\n"
"    for (cc=0; cc<ncol_z; ++cc) {\n"
"      /* Get the dense column of y */\n"
"      for (kk=colind_y[cc]; kk<colind_y[cc+1]; ++kk) {\n"
"        w[row_y[kk]] = y[kk];\n"
"      }\n"
"      /* Loop over the nonzeros of z */\n"
"      for (kk=colind_z[cc]; kk<colind_z[cc+1]; ++kk) {\n"
"        int rr = row_z[kk];\n"
"        /* Loop over corresponding columns of x */\n"
"        for (kk1=colind_x[rr]; kk1<colind_x[rr+1]; ++kk1) {\n"
"          z[kk] += x[kk1] * w[row_x[kk1]];\n"
"        }\n"
"      }\n"
"    }\n"
"  } else {\n"
"    /* Loop over the columns of y and z */\n"
"    for (cc=0; cc<ncol_y; ++cc) {\n"
"      /* Get the dense column of z */\n"
"      for (kk=colind_z[cc]; kk<colind_z[cc+1]; ++kk) {\n"
"        w[row_z[kk]] = z[kk];\n"
"      }\n"
"      /* Loop over the nonzeros of y */\n"
"      for (kk=colind_y[cc]; kk<colind_y[cc+1]; ++kk) {\n"
"        int rr = row_y[kk];\n"
"        /* Loop over corresponding columns of x */\n"
"        for (kk1=colind_x[rr]; kk1<colind_x[rr+1]; ++kk1) {\n"
"          w[row_x[kk1]] += x[kk1]*y[kk];\n"
"        }\n"
"      }\n"
"      /* Get the sparse column of z */\n"
"      for (kk=colind_z[cc]; kk<colind_z[cc+1]; ++kk) {\n"
"        z[kk] = w[row_z[kk]];\n"
"      }\n"
"    }\n"
"  }\n"
"}\n";

const char * codegen_str_mv =
"void CASADI_PREFIX(mv)(const real_t* x, const int* sp_x, const real_t* y, real_t* z, int tr) {\n"
"  /* Get sparsities */\n"
"  int ncol_x = sp_x[1];\n"
"  const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;\n"
"  int i, el;\n"
"  if (tr) {\n"
"    /* loop over the columns of x */\n"
"    for (i=0; i<ncol_x; ++i) {\n"
"      /* loop over the non-zeros of x */\n"
"      for (el=colind_x[i]; el<colind_x[i+1]; ++el) {\n"
"        z[i] += x[el] * y[row_x[el]];\n"
"      }\n"
"    }\n"
"  } else {\n"
"    /* loop over the columns of x */\n"
"    for (i=0; i<ncol_x; ++i) {\n"
"      /* loop over the non-zeros of x */\n"
"      for (el=colind_x[i]; el<colind_x[i+1]; ++el) {\n"
"        z[row_x[el]] += x[el] * y[i];\n"
"      }\n"
"    }\n"
"  }\n"
"}\n";

const char * codegen_str_trans =
"void CASADI_PREFIX(trans)(const real_t* x, const int* sp_x, real_t* y, const int* sp_y, int *tmp) {\n"
"  int ncol_x = sp_x[1];\n"
"  int nnz_x = sp_x[2 + ncol_x];\n"
"  const int* row_x = sp_x + 2 + ncol_x+1;\n"
"  int ncol_y = sp_y[1];\n"
"  const int* colind_y = sp_y+2;\n"
"  int k;\n"
"  for (k=0; k<ncol_y; ++k) tmp[k] = colind_y[k];\n"
"  for (k=0; k<nnz_x; ++k) {\n"
"    y[tmp[row_x[k]]++] = x[k];\n"
"  }\n"
"}\n";

const char * codegen_str_norm_1 =
"real_t CASADI_PREFIX(norm_1)(int n, const real_t* x) {\n"
"  real_t ret = 0;\n"
"  int i;\n"
"  if (x) {\n"
"    for (i=0; i<n; ++i) ret += fabs(*x++);\n"
"  }\n"
"  return ret;\n"
"}\n";

const char * codegen_str_norm_2 =
"real_t CASADI_PREFIX(norm_2)(int n, const real_t* x) {\n"
"  return sqrt(CASADI_PREFIX(dot)(n, x, x));\n"
"}\n";

const char * codegen_str_norm_inf =
"real_t CASADI_PREFIX(norm_inf)(int n, const real_t* x) {\n"
"  real_t ret = 0;\n"
"  int i;\n"
"  for (i=0; i<n; ++i) ret = fmax(ret, fabs(*x++));\n"
"  return ret;\n"
"}\n";

const char * codegen_str_norm_inf_mul =
"real_t CASADI_PREFIX(norm_inf_mul)(const real_t* x, const int* sp_x, const real_t* y, const int* sp_y, real_t* dwork, int* iwork) {\n"
"  real_t res = 0;\n"
"  /* Get sparsities */\n"
"  int nrow_x = sp_x[0], ncol_x = sp_x[1];\n"
"  const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;\n"
"  int ncol_y = sp_y[1];\n"
"  const int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;\n"
"  /* Implementation borrowed from Scipy's sparsetools/csr.h */\n"
"  /* method that uses O(n) temp storage */\n"
"  int *mask = iwork + ncol_y+1;\n"
"  int i,jj,kk;\n"
"  // Pass 1\n"
"  for (i=0; i<nrow_x; ++i) mask[i] = -1;\n"
"  iwork[0] = 0;\n"
"  int nnz = 0;\n"
"  for (i=0; i<ncol_y; ++i) {\n"
"    int row_nnz = 0;\n"
"    for (jj=colind_y[i]; jj < colind_y[i+1]; jj++) {\n"
"      int j = row_y[jj];\n"
"      for (kk=colind_x[j]; kk < colind_x[j+1]; kk++) {\n"
"        int k = row_x[kk];\n"
"        if (mask[k] != i) {\n"
"          mask[k] = i;\n"
"          row_nnz++;\n"
"        }\n"
"      }\n"
"    }\n"
"    int next_nnz = nnz + row_nnz;\n"
"    nnz = next_nnz;\n"
"    iwork[i+1] = nnz;\n"
"  }\n"
"  // Pass 2\n"
"  int *next = iwork + ncol_y+1;\n"
"  for (i=0; i<nrow_x; ++i) next[i] = -1;\n"
"  real_t* sums = dwork;\n"
"  for (i=0; i<nrow_x; ++i) sums[i] = 0;\n"
"  nnz = 0;\n"
"  iwork[0] = 0;\n"
"  for (i=0; i<ncol_y; ++i) {\n"
"    int head   = -2;\n"
"    int length =  0;\n"
"    int jj_start = colind_y[i];\n"
"    int jj_end   = colind_y[i+1];\n"
"    for (jj=jj_start; jj<jj_end; ++jj) {\n"
"      int j = row_y[jj];\n"
"      real_t v = y[jj];\n"
"      int kk_start = colind_x[j];\n"
"      int kk_end   = colind_x[j+1];\n"
"      for (kk = kk_start; kk<kk_end; ++kk) {\n"
"        int k = row_x[kk];\n"
"        sums[k] += v*x[kk];\n"
"        if (next[k] == -1) {\n"
"          next[k] = head;\n"
"          head  = k;\n"
"          length++;\n"
"        }\n"
"      }\n"
"    }\n"
"    for (jj=0; jj<length; ++jj) {\n"
"      if (!is_zero(sums[head])) {\n"
"        real_t a = fabs(sums[head]);\n"
"        res = fmax(res, a);\n"
"        nnz++;\n"
"      }\n"
"      int temp = head;\n"
"      head = next[head];\n"
"      next[temp] = -1; //clear arrays\n"
"      sums[temp] =  0;\n"
"    }\n"
"    iwork[i+1] = nnz;\n"
"  }\n"
"  return res;\n"
"}\n";

const char * codegen_str_bilin =
"real_t CASADI_PREFIX(bilin)(const real_t* A, const int* sp_A, const real_t* x, const real_t* y) {\n"
"  /* Get sparsities */\n"
"  int ncol_A = sp_A[1];\n"
"  const int *colind_A = sp_A+2, *row_A = sp_A + 2 + ncol_A+1;\n"
"  /* Return value */\n"
"  real_t ret=0;\n"
"  /* Loop over the columns of A */\n"
"  int cc, rr, el;\n"
"  for (cc=0; cc<ncol_A; ++cc) {\n"
"    /* Loop over the nonzeros of A */\n"
"    for (el=colind_A[cc]; el<colind_A[cc+1]; ++el) {\n"
"      /* Get the row */\n"
"      rr=row_A[el];\n"
"      /* Add contribution */\n"
"      ret += x[rr]*A[el]*y[cc];\n"
"    }\n"
"  }\n"
"  return ret;\n"
"}\n";

const char * codegen_str_rank1 =
"void CASADI_PREFIX(rank1)(real_t* A, const int* sp_A, real_t alpha, const real_t* x, const real_t* y) {\n"
"  /* Get sparsities */\n"
"  int ncol_A = sp_A[1];\n"
"  const int *colind_A = sp_A+2, *row_A = sp_A + 2 + ncol_A+1;\n"
"  /* Loop over the columns of A */\n"
"  int cc, rr, el;\n"
"  for (cc=0; cc<ncol_A; ++cc) {\n"
"    /* Loop over the nonzeros of A */\n"
"    for (el=colind_A[cc]; el<colind_A[cc+1]; ++el) {\n"
"      /* Get row */\n"
"      rr = row_A[el];\n"
"      /* Add the multiple */\n"
"      A[el] += alpha*x[rr]*y[cc];\n"
"    }\n"
"  }\n"
"}\n";

const char * codegen_str_getu =
"void CASADI_PREFIX(getu)(const real_t* x, const int* sp_x, real_t* v) {\n"
"  /* Get sparsities */\n"
"  int ncol_x = sp_x[1];\n"
"  const int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;\n"
"  /* Loop over the columns of x */\n"
"  int cc, el;\n"
"  for (cc=0; cc<ncol_x; ++cc) {\n"
"    /* Loop over the nonzeros of x */\n"
"    for (el=colind_x[cc]; el<colind_x[cc+1] && row_x[el]<=cc; ++el) {\n"
"      *v++ = x[el];\n"
"    }\n"
"  }\n"
"}\n";

const char * codegen_str_polyval =
"real_t CASADI_PREFIX(polyval)(const real_t* p, int n, real_t x) {\n"
"  real_t r=p[0];\n"
"  int i;\n"
"  for (i=1; i<=n; ++i) {\n"
"    r = r*x + p[i];\n"
"  }\n"
"  return r;\n"
"}\n";

const char * codegen_str_flip =
"int CASADI_PREFIX(flip)(int* corner, int ndim) {\n"
"  int i;\n"
"  for (i=0; i<ndim; ++i) {\n"
"    if (corner[i]) {\n"
"      corner[i]=0;\n"
"    } else {\n"
"      corner[i]=1;\n"
"      return 1;\n"
"    }\n"
"  }\n"
"  return 0;\n"
"}\n";

const char * codegen_str_low =
"int CASADI_PREFIX(low)(real_t x, const double* grid, int ng) {\n"
"  int i;\n"
"  for (i=0; i<ng-2; ++i) {\n"
"    if (x < grid[i+1]) break;\n"
"  }\n"
"  return i;\n"
"}\n";

const char * codegen_str_interpn_weights =
"void CASADI_PREFIX(interpn_weights)(int ndim, const real_t* grid, const int* offset, const real_t* x, real_t* alpha, int* index) {\n"
"  /* Left index and fraction of interval */\n"
"  int i;\n"
"  for (i=0; i<ndim; ++i) {\n"
"    /* Grid point */\n"
"    real_t xi = x ? x[i] : 0;\n"
"    /* Grid */\n"
"    const real_t* g = grid + offset[i];\n"
"    int ng = offset[i+1]-offset[i];\n"
"    /* Find left index */\n"
"    int j = index[i] = CASADI_PREFIX(low)(xi, g, ng);\n"
"    /* Get interpolation/extrapolation alpha */\n"
"    alpha[i] = (xi-g[j])/(g[j+1]-g[j]);\n"
"  }\n"
"}\n";

const char * codegen_str_interpn_interpolate =
"real_t CASADI_PREFIX(interpn_interpolate)(int ndim, const int* offset, const real_t* values, const real_t* alpha, const int* index, const int* corner, real_t* coeff) {\n"
"  /* Get weight and value for corner */\n"
"  real_t c=1;\n"
"  int ld=1; /* leading dimension */\n"
"  int i;\n"
"  for (i=0; i<ndim; ++i) {\n"
"    if (coeff) *coeff++ = c;\n"
"    if (corner[i]) {\n"
"      c *= alpha[i];\n"
"    } else {\n"
"      c *= 1-alpha[i];\n"
"    }\n"
"    values += (index[i]+corner[i])*ld;\n"
"    ld *= offset[i+1]-offset[i];\n"
"  }\n"
"  if (coeff) {\n"
"    return *values;\n"
"  } else {\n"
"    return c**values;\n"
"  }\n"
"}\n";

const char * codegen_str_interpn =
"real_t CASADI_PREFIX(interpn)(int ndim, const real_t* grid, const int* offset, const real_t* values, const real_t* x, int* iw, real_t* w) {\n"
"  /* Work vectors */\n"
"  real_t* alpha = w; w += ndim;\n"
"  int* index = iw; iw += ndim;\n"
"  int* corner = iw; iw += ndim;\n"
"  /* Left index and fraction of interval */\n"
"  CASADI_PREFIX(interpn_weights)(ndim, grid, offset, x, alpha, index);\n"
"  /* Loop over all corners, add contribution to output */\n"
"  CASADI_PREFIX(fill_int)(corner, ndim, 0);\n"
"  real_t ret = 0;\n"
"  do {\n"
"    real_t* coeff = 0;\n"
"    ret += CASADI_PREFIX(interpn_interpolate)(ndim, offset, values,\n"
"      alpha, index, corner, coeff);\n"
"  } while (CASADI_PREFIX(flip)(corner, ndim));\n"
"  return ret;\n"
"}\n";

const char * codegen_str_interpn_grad =
"void CASADI_PREFIX(interpn_grad)(real_t* grad, int ndim, const real_t* grid, const int* offset, const real_t* values, const real_t* x, int* iw, real_t* w) {\n"
"  /* Quick return */\n"
"  if (!grad) return;\n"
"  /* Work vectors */\n"
"  real_t* alpha = w; w += ndim;\n"
"  real_t* coeff = w; w += ndim;\n"
"  int* index = iw; iw += ndim;\n"
"  int* corner = iw; iw += ndim;\n"
"  /* Left index and fraction of interval */\n"
"  CASADI_PREFIX(interpn_weights)(ndim, grid, offset, x, alpha, index);\n"
"  /* Loop over all corners, add contribution to output */\n"
"  CASADI_PREFIX(fill_int)(corner, ndim, 0);\n"
"  CASADI_PREFIX(fill)(grad, ndim, 0.);\n"
"  do {\n"
"    /* Get coefficients */\n"
"    real_t v = CASADI_PREFIX(interpn_interpolate)(ndim, offset, values,\n"
"      alpha, index, corner, coeff);\n"
"    /* Propagate to alpha */\n"
"    int i;\n"
"    for (i=ndim-1; i>=0; --i) {\n"
"      if (corner[i]) {\n"
"        grad[i] += v*coeff[i];\n"
"        v *= alpha[i];\n"
"      } else {\n"
"        grad[i] -= v*coeff[i];\n"
"        v *= 1-alpha[i];\n"
"      }\n"
"    }\n"
"  } while (CASADI_PREFIX(flip)(corner, ndim));\n"
"  /* Propagate to x */\n"
"  int i;\n"
"  for (i=0; i<ndim; ++i) {\n"
"    const real_t* g = grid + offset[i];\n"
"    int j = index[i];\n"
"    grad[i] /= g[j+1]-g[j];\n"
"  }\n"
"}\n";

const char * codegen_str_copy_define = "#define copy(x, n, y) CASADI_PREFIX(copy)(x, n, y)\n";
const char * codegen_str_swap_define = "#define swap(n, x, inc_x, y, inc_y) CASADI_PREFIX(swap)(n, x, inc_x, y, inc_y)\n";
const char * codegen_str_project_define = "#define project(x, sp_x, y, sp_y, w) CASADI_PREFIX(project)(x, sp_x, y, sp_y, w)\n";
const char * codegen_str_mproject_define = "#define mproject(factor, x, sp_x, y, sp_y, w) CASADI_PREFIX(mproject)(factor, x, sp_x, y, sp_y, w)\n";
const char * codegen_str_maddproject_define = "#define maddproject(factor, x, sp_x, y, sp_y, w) CASADI_PREFIX(maddproject)(factor, x, sp_x, y, sp_y, w)\n";
const char * codegen_str_dense_transfer_define = "#define dense_transfer(factor, x, sp_x, y, sp_y, w) CASADI_PREFIX(dense_transfer)(factor, x, sp_x, y, sp_y, w)\n";
const char * codegen_str_densify_define = "#define densify(x, sp_x, y, tr) CASADI_PREFIX(densify)(x, sp_x, y, tr)\n";
const char * codegen_str_sparsify_define = "#define sparsify(x, y, sp_y, tr) CASADI_PREFIX(sparsify)(x, y, sp_y, tr)\n";
const char * codegen_str_scal_define = "#define scal(n, alpha, x) CASADI_PREFIX(scal)(n, alpha, x)\n";
const char * codegen_str_axpy_define = "#define axpy(n, alpha, x, y) CASADI_PREFIX(axpy)(n, alpha, x, y)\n";
const char * codegen_str_dot_define = "#define dot(n, x, y) CASADI_PREFIX(dot)(n, x, y)\n";
const char * codegen_str_max_viol_define = "#define max_viol(n, x, lb, ub) CASADI_PREFIX(max_viol)(n, x, lb, ub)\n";
const char * codegen_str_sum_viol_define = "#define sum_viol(n, x, lb, ub) CASADI_PREFIX(sum_viol)(n, x, lb, ub)\n";
const char * codegen_str_iamax_define = "#define iamax(n, x, inc_x) CASADI_PREFIX(iamax)(n, x, inc_x)\n";
const char * codegen_str_fill_define = "#define fill(x, n, alpha) CASADI_PREFIX(fill)(x, n, alpha)\n";
const char * codegen_str_fill_int_define = "#define fill_int(x, n, alpha) CASADI_PREFIX(fill_int)(x, n, alpha)\n";
const char * codegen_str_mtimes_define = "#define mtimes(x, sp_x, y, sp_y, z, sp_z, w, tr) CASADI_PREFIX(mtimes)(x, sp_x, y, sp_y, z, sp_z, w, tr)\n";
const char * codegen_str_mv_define = "#define mv(x, sp_x, y, z, tr) CASADI_PREFIX(mv)(x, sp_x, y, z, tr)\n";
const char * codegen_str_trans_define = "#define trans(x, sp_x, y, sp_y,  int *tmp) CASADI_PREFIX(trans)(x, sp_x, y, sp_y,  int *tmp)\n";
const char * codegen_str_norm_1_define = "#define norm_1(n, x) CASADI_PREFIX(norm_1)(n, x)\n";
const char * codegen_str_norm_2_define = "#define norm_2(n, x) CASADI_PREFIX(norm_2)(n, x)\n";
const char * codegen_str_norm_inf_define = "#define norm_inf(n, x) CASADI_PREFIX(norm_inf)(n, x)\n";
const char * codegen_str_norm_inf_mul_define = "#define norm_inf_mul(x, sp_x, y, sp_y, dwork, iwork) CASADI_PREFIX(norm_inf_mul)(x, sp_x, y, sp_y, dwork, iwork)\n";
const char * codegen_str_bilin_define = "#define bilin(A, sp_A, x, y) CASADI_PREFIX(bilin)(A, sp_A, x, y)\n";
const char * codegen_str_rank1_define = "#define rank1(A, sp_A, alpha, x, y) CASADI_PREFIX(rank1)(A, sp_A, alpha, x, y)\n";
const char * codegen_str_getu_define = "#define getu(x, sp_x, v) CASADI_PREFIX(getu)(x, sp_x, v)\n";
const char * codegen_str_polyval_define = "#define polyval(p, n, x) CASADI_PREFIX(polyval)(p, n, x)\n";
const char * codegen_str_flip_define = "#define flip(corner, ndim) CASADI_PREFIX(flip)(corner, ndim)\n";
const char * codegen_str_low_define = "#define low(x, grid, ng) CASADI_PREFIX(low)(x, grid, ng)\n";
const char * codegen_str_interpn_weights_define = "#define interpn_weights(ndim, grid, offset, x, alpha, index) CASADI_PREFIX(interpn_weights)(ndim, grid, offset, x, alpha, index)\n";
const char * codegen_str_interpn_interpolate_define = "#define interpn_interpolate(ndim, offset, values, alpha, index, corner, coeff) CASADI_PREFIX(interpn_interpolate)(ndim, offset, values, alpha, index, corner, coeff)\n";
const char * codegen_str_interpn_define = "#define interpn(ndim, grid, offset, values, x, iw, w) CASADI_PREFIX(interpn)(ndim, grid, offset, values, x, iw, w)\n";
const char * codegen_str_interpn_grad_define = "#define interpn_grad(grad, ndim, grid, offset, values, x, iw, w) CASADI_PREFIX(interpn_grad)(grad, ndim, grid, offset, values, x, iw, w)\n";
}