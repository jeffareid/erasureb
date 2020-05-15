//----------------------------------------------------------------------//
//      ersbch.cpp  erasure demo - BCH based RS code                    //
//                  Copyright(c) 2020, Jeff Reid                        //
//                  2020MAY15 14:30                                     //
//----------------------------------------------------------------------//
//      equates                                                         //
//----------------------------------------------------------------------//
#include <ctime>
#include <iostream>
#include <iomanip>
#include <memory>

typedef unsigned char      BYTE;
typedef unsigned long long QWORD;

#define NDAT    17                      // # data rows
#define NPAR    3                       // # parity rows
#define NCOL    (32*1024)               // # columns (multiple of 16)

#define DISPLAYI 0                      // dislay matrixinv

const int NROW = NDAT+NPAR;             // # rows (total)

class MATRIX{
public:
    BYTE *p = NULL;                     // ptr to data
    BYTE *m[NROW];                      // array of ptrs
    int r;                              // # rows
    int c;                              // # columns

MATRIX(int rr, int cc)                  // new matrix
  {
    p = new BYTE[rr*cc];
    r = rr;
    c = cc;
    for(int i = 0; i < rr; i++)
        m[i] = p+i*cc;
  }

MATRIX(int rr, int cc, BYTE *pp)        // mapped matrix
  {
    p = NULL;
    r = rr;
    c = cc;
    for(int i = 0; i < rr; i++)
        m[i] = pp+i*cc;
  }

~MATRIX()                               // destructor
  {
    if(p != NULL)
        delete[] p;
  }  
};

//----------------------------------------------------------------------//
//      data                                                            //
//----------------------------------------------------------------------//
static int aiGA[60] =                   // GF's and Alpha's //
   {0x11b,0x03,0x11d,0x02,0x12b,0x02,0x12d,0x02,
    0x139,0x03,0x13f,0x03,0x14d,0x02,0x15f,0x02,
    0x163,0x02,0x165,0x02,0x169,0x02,0x171,0x02,
    0x177,0x03,0x17b,0x0b,0x187,0x02,0x18b,0x0b,
    0x18d,0x02,0x19f,0x03,0x1a3,0x03,0x1a9,0x02,
    0x1b1,0x07,0x1bd,0x07,0x1c3,0x02,0x1cf,0x02,
    0x1d7,0x07,0x1dd,0x07,0x1e7,0x02,0x1f3,0x0d,
    0x1f5,0x02,0x1f9,0x03};

static BYTE     abExp[512];             // antilog table
static BYTE     abLog[256];             // log table

static BYTE     abId[512];              // used by MatrixInv

static BYTE     abMpy[65536];           // matrix for GFMpy
static BYTE     abDiv[65536];           // matrix for GFDiv
static BYTE     abPow[65536];           // matrix for GFPow

static BYTE     abRoot[NPAR];
static BYTE     abPoly[NPAR+1];

static int      iGF;                    // Galios Field Polynomial
static BYTE     bAlpha;                 // Alpha for this field

static clock_t  ctTimeStart;            // clock values
static clock_t  ctTimeStop;

//----------------------------------------------------------------------//
//      code                                                            //
//----------------------------------------------------------------------//
#define GFAdd(a, b) ((a)^(b))
#define GFSub(a, b) ((a)^(b))
#define GFMpy(a, b) (abMpy[((a)<<8)+(b)])
#define GFDiv(a, b) (abDiv[((a)<<8)+(b)])
#define GFPow(a, b) (abPow[((a)<<8)+(b)])

static BYTE     GFMpy0(BYTE, BYTE);
static BYTE     GFMpy1(BYTE, BYTE);
static BYTE     GFDiv1(BYTE, BYTE);
static BYTE     GFPow1(BYTE, BYTE);
static void     InitGF(void);

static void     MatrixMpy(MATRIX &, MATRIX &, MATRIX &);
static int      MatrixInv(MATRIX &, MATRIX &);
static void     MatrixXor(MATRIX &, int);
static void     ShowMatrix(MATRIX &);

static void     InitCombination(int[], int k, int n);
static bool     NextCombination(int [], int, int);
static void     Patterns(MATRIX &, MATRIX &);

//----------------------------------------------------------------------//
//      GFMpy0(b0,b1)           b0*b1       using low level math        //
//----------------------------------------------------------------------//
static BYTE GFMpy0(BYTE b0, BYTE b1)
{
int i;
int product;
    product = 0;
    for(i = 0; i < 8; i++){
        product <<= 1;
        if(product & 0x100){
            product ^= iGF;}
        if(b0 & 0x80u){
            product ^= b1;}
        b0 <<= 1;}
    return((BYTE)product);
}

//----------------------------------------------------------------------//
//      GFMpy1(b0, b1)           b0*b1       using logs                 //
//----------------------------------------------------------------------//
static BYTE GFMpy1(BYTE byt0, BYTE byt1)
{
    if(byt0 == 0 || byt1 == 0)
        return(0);

    return(abExp[(int)abLog[byt0]+(int)abLog[byt1]]);
}

//----------------------------------------------------------------------//
//      GFDiv1(b0, b1)           b0/b1                                  //
//----------------------------------------------------------------------//
static BYTE GFDiv1(BYTE b0, BYTE b1)
{
    if(b1 == 0){
        std::cout << "divide by zero" << std::endl;
        return(0);}
    if(b0 == 0)
        return(0);
    return(abExp[255+(int)abLog[b0]-(int)abLog[b1]]);
}

//----------------------------------------------------------------------//
//      GFPow1(b0, b1)           b0^b1                                  //
//----------------------------------------------------------------------//
static BYTE GFPow1(BYTE b0, BYTE b1)
{
BYTE b;
    b = 1;
    while(b1){
        if(b1&1)
            b = GFMpy(b, b0);
        b0 = GFMpy(b0, b0);
        b1 >>= 1;}
    return(b);
}

//----------------------------------------------------------------------//
//      InitGF  Initialize Galios Stuff                                 //
//----------------------------------------------------------------------//
static void InitGF(void)
{
BYTE b;
int i, j;

    b = 1;
    for(i = 0; i < 512; i++){           // init abExp[] //
        abExp[i] = b;
        b = GFMpy0(b, bAlpha);}

    abLog[0] = 0xff;                    // init abLog[] //
    for(i = 0; i < 255; i++){
        abLog[abExp[i]] = i;}

    for(j = 0; j < 256; j++){           // init GF math matrices
        for (i = 0; i < 256; i++){
            abMpy[(j<<8)+i] = GFMpy1(j, i);
            if(i != 0){
                abDiv[(j<<8)+i] = GFDiv1(j, i);}
        }
    }
    for (j = 0; j < 256; j++) {
        for (i = 0; i < 256; i++) {
            abPow[(j << 8) + i] = GFPow1(j, i);
        }
    }

    memset(abId, 0, sizeof(abId));      // set up abId for MatrixInv
    abId[256] = 1;

    for(i = 0; i < NPAR; i++)           // init abRoot
        abRoot[i] = GFPow(2, i);

    abPoly[0] = 1;                      // init abPoly
    for(j = 0; j < NPAR; j++){
        for(i = j; i >= 0; i--){
            abPoly[i+1] = GFSub(abPoly[i+1],
                    GFMpy(abPoly[i], abRoot[j]));}}
}

//----------------------------------------------------------------------//
//      MatrixMpy(mDst, mSrc0, mSrc1) matrix multiply                   //
//----------------------------------------------------------------------//
static void MatrixMpy(MATRIX &mDst, MATRIX &mSrc0, MATRIX &mSrc1)
{
int r, m, c;

    for(r = 0; r < mSrc0.r; r++){       // for each row
        memset(mDst.m[r], 0, mDst.c);   // zero dst row
        for(m = 0; m < mSrc0.c; m++){   // do the mpy
            for (c = 0; c < mSrc1.c; c++){
                mDst.m[r][c] ^= GFMpy(mSrc0.m[r][m], mSrc1.m[m][c]);
            }
        }
    }
}

//----------------------------------------------------------------------//
//      MatrixInv(mDst, mSrc) invert matrix                             //
//      assumes square matrix                                           //
//----------------------------------------------------------------------//
static int MatrixInv(MATRIX &mDst, MATRIX &mSrc)
{
int r, c, t;
BYTE b;

    MATRIX mAug(mSrc.r, mSrc.c<<1);     // create augmented matrix

//      generate augmented matrix
//      left side is copy of mSrc, right side is Identity Matrix

    for(r = 0; r < mSrc.r; r++){
        memcpy(mAug.m[r],        mSrc.m[r],  mSrc.c);   // copy row of data
        memcpy(mAug.m[r]+mSrc.c, abId+256-r, mSrc.c);   // append id
    }

//      normalize according to left side
//      results in inverse matrix in right side

    #if DISPLAYI
        printf("start\n");
        ShowMatrix(mAug);
    #endif

    for(r = 0; r < mSrc.r; r++){            // working at [r][r]
//      find 1st non-zero in current row's column
//      (no check for non-invertible matrix)
        t = r;
        while(0 == mAug.m[t][r])
            t++;
//      swap rows if needed
        if(t != r){
            std::swap(mAug.m[r], mAug.m[t]);
            #if DISPLAYI
                printf("swapped rows\n");
                ShowMatrix(mAug);
            #endif
        }
//      divide row by [r][r] so [r][r] == 1
        b = mAug.m[r][r];
        for(c = r; c < mAug.c; c++)
            mAug.m[r][c] = GFDiv(mAug.m[r][c], b);
        #if DISPLAYI
            printf("divided row\n");
            ShowMatrix(mAug);
        #endif
//      subtract multiple of this row from other rows
//      to create a column of zeroes
        for(t = 0; t < mAug.r; t++){
            if(t == r)                      // skip if current row
                continue;
            b = mAug.m[t][r];
            for(c = r; c < mAug.c; c++)
                mAug.m[t][c] = GFSub(mAug.m[t][c], GFMpy(b, mAug.m[r][c]));
        }
#if DISPLAYI
        printf("zeroed columns\n");
        ShowMatrix(mAug);
#endif
    }
//      copy right side of mAug to mDst
    for(r = 0; r < mDst.r; r++)
        memcpy(mDst.m[r], mAug.m[r]+mSrc.c, mSrc.c);
    return(0);
}

//----------------------------------------------------------------------//
//      MatrixXor - xor all but [row][] to [row]                        //
//----------------------------------------------------------------------//
static void MatrixXor(MATRIX &mDst, int row)
{
int r, c;
BYTE *pSrc;                                 // ptr to src row
    pSrc = mDst.m[0]; 
    r = 1;
    if(row == 0){
        pSrc = mDst.m[1]; 
        r = 2;
    }
    for(c = 0; c < NCOL; c += 8)            // copy src row to dst row
        *(QWORD *)&mDst.m[row][c] = *(QWORD *)(pSrc+c);
    for( ; r < NROW; r++){                  // xor  non-erased rows
        if(r == row)
            continue;
        for(c = 0; c < NCOL; c += 8)
            *(QWORD *)&mDst.m[row][c] ^= *(QWORD *)&mDst.m[r][c];
    }
}

//----------------------------------------------------------------------//
//      ShowMatrix                                                      //
//----------------------------------------------------------------------//
static void ShowMatrix(MATRIX &mSrc)
{
int  r, c;
    for(r = 0; r < mSrc.r; r++){
        for(c = 0; c < mSrc.c; c++)
            std::cout << std::setfill('0') << std::hex
                      << std::setw(2) << (int)mSrc.m[r][c];
        std::cout << std::endl;}
    std::cout << std::endl;
}

//----------------------------------------------------------------------//
//      InitCombination - init combination to first set - 1             //
//----------------------------------------------------------------------//
void InitCombination(int a[], int k, int n) {
    for(int i = 0; i < k; i++)
        a[i] = i;
    --a[k-1];     // 1st call to NextCombination will return 1st set
}

//----------------------------------------------------------------------//
//      NextCombination - generate next combination                     //
//----------------------------------------------------------------------//
bool NextCombination(int a[], int k, int n) {
int j = k - 1;
    while (j >= 0 && a[j] == n - k + j)
        --j;
    if (j == -1)
        return false;
    ++a[j];
    for (int i = j + 1; i < k; ++i)
        a[i] = a[j] + i - j;
    return true;
}

//----------------------------------------------------------------------//
//      Patterns   test all erasure patterns                            //
//----------------------------------------------------------------------//
static void Patterns(MATRIX &mSyn, MATRIX &mDat)
{
int e[NPAR];                                // erasures
int r, c, n, m;

    for(r = 0; r < NROW; r++){              // single erasures
        memset(mDat.m[r], 0xaa, NCOL);      // corrupt row
        MatrixXor(mDat, r);                 // xor non-erased rows
    }
        
    for(n = 2; n <= NPAR; n++){             // n = number of erasures
        m = n - 1;                          // m = n-1
        MATRIX mLct(n,n);                   // locator matrix
        MATRIX mInv(n,n);                   // inverse matrix
        MATRIX mCor(m,NROW);                // correction matrix
        MATRIX mFix(m,NCOL,mDat.m[0]);      // fixed data matrix
        InitCombination(e, n, NROW);        // setup next combination
        while(NextCombination(e, n, NROW)){ // set e == erasure indexes
            for(r = 0; r < n; r++)          // set mFix to mDat erased rows
                mFix.m[r] = mDat.m[e[r]];
            for(r = 0; r < n; r++)          // corrupt erased rows
                memset(mDat.m[e[r]], 0xaa, NCOL);
            for(r = 0; r < n; r++)          // generate locator matrix
                for(c = 0; c < n; c++)
                    mLct.m[r][c] = GFPow(abRoot[r], NROW-1-e[c]);
            MatrixInv(mInv, mLct);          // invert locator matrix
            mInv.r = m;                     // change mInv to m rows
            MatrixMpy(mCor, mInv, mSyn);    // generate m correction rows
            for(r = 0; r < m; r++)          // zero out n erased columns
                for(c = 0; c < n; c++)
                    mCor.m[r][e[c]] = 0;
            MatrixMpy(mFix, mCor, mDat);    // correct m rows
            MatrixXor(mDat, e[m]);          // correct last row via xor
            mInv.r = n;                     // restore mInv to n rows
        }
    }
}

//----------------------------------------------------------------------//
//      main                                                            //
//----------------------------------------------------------------------//
int main()
{
int r, c;
BYTE b;

    iGF = aiGA[2];                          // select GF params
    bAlpha = aiGA[3];
    InitGF();
    MATRIX mEnc(NPAR, NDAT);                // generate encode matrix
    {
        BYTE abRem[NPAR+1] = {0};
        BYTE q = 1;
        for(c = NDAT-1; c >= 0; c--){
            for(int i = 0; i < NPAR; i++)
                abRem[i] = GFSub(abRem[i+1], GFMpy(q, abPoly[i+1]));
            for(r = 0; r < NPAR; r++)
                mEnc.m[r][c] = abRem[r];
            q = abRem[0];
        }
    }

    MATRIX mSyn(NPAR, NROW);               // generate syndrome matrix
    for(r = 0; r < NPAR; r++)
        for(c = 0; c < NROW; c++)
            mSyn.m[r][c] = GFPow(2,r*(NROW-1-c));

    MATRIX mDat(NROW, NCOL);               // generate data matrix
    b = 0;
    for (r = 0; r < NDAT; r++)
        for (c = 0; c < NCOL; c++)
            mDat.m[r][c] = b++;

    // generate parity matrix mapped into data matrix
    MATRIX mPar(NPAR, NCOL, mDat.m[0]+NDAT*NCOL);
    MatrixMpy(mPar, mEnc, mDat);           // encode data

    ctTimeStart = clock();
    Patterns(mSyn, mDat);                  // test erasure patterns
    ctTimeStop = clock();
    std::cout << "# of ticks " << ctTimeStop - ctTimeStart << std::endl;

    b = 0;                                 // do a one time verify of mDat
    for (r = 0; r < NDAT; r++) {
        for (c = 0; c < NCOL; c++) {
            if (mDat.m[r][c] != b++) {
                goto vfy0;}}}
vfy0:
    if(r == NDAT)
        std::cout << "passed" << std::endl;
    else
        std::cout << "failed" << std::endl;

    return 0;
}
