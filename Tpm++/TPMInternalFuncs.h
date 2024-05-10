#pragma once
void fmgN(CDGMulDimData &A, CDGMulDimData &b, double p[8], double* p4, CDGMulDimData& rout);
void mrf1(unsigned int dm[], unsigned char q[], float p[], float G[], float w[], int code);
bool AnyWarpTissueSetting(sTissues & tissues, bool NativeTissue_NativeSpace, bool NativeTissue_DartelImported, bool WarpTissue_Modulated, bool WarpTissue_UnModulated);
void ZeroImageMem(sitk::Image & img);
MAT_TYPE BSplines(double *c, DimVec dim, MAT_TYPE &x0, MAT_TYPE &x1, MAT_TYPE &x2, DoubleVec &deg);
bool SaveStringToFile(string s, string filepath);
string LoadStringFromFile(string filepath);
const double Div255d = 1.0 / 255;
void Clean_Gwc(CDGMulDimData &  PBytes, int level);
int Slice1(double * mat, float* image, int xdim1, int ydim1, CDGMulDimData& vol, double background, DoubleVec & scale, DoubleVec &offset);
VECTOR_TYPE VoxelSizeSigned(MAT_TYPE & mat);
void GetBoundBox(MAT_TYPE & tpmmat, SizeVec &tpmDim, MAT_TYPE & bb, VECTOR_TYPE &vx);
void MarkrofField(IN OUT CDGMulDimData & vPByte, IN CDGMulDimData& QFl, IN CDGMulDimData& vG, IN VECTOR_TYPE &vx2);
void RGrid(DimVec d, CDGMulDimData & x);
void affind(CDGMulDimData & y0, CDGMulDimData & y1, MAT_TYPE & M);
void GetClosestAffine(CDGMulDimData & x, CDGMulDimData & y, CDGMulDimData &w1, CDGMulDimData &w2, OUT MAT_TYPE &M, OUT MAT_TYPE *R);
void ExtrapolateDef(CDGMulDimData &Y, MAT_TYPE & M);
void defsw(IN vector<CDGMulDimData>& Coef, IN  int z, IN MAT_TYPE &x0, IN MAT_TYPE &y0, IN VECTOR_TYPE &z0, IN MAT_TYPE &M, IN MAT_TYPE &MT, DoubleVec & prm,
	OUT MAT_TYPE &x1, OUT MAT_TYPE &y1, OUT MAT_TYPE &z1);
void likelihoods(vector<MAT_TYPE> &f, vector<MAT_TYPE> &bf, VECTOR_TYPE& mg, MAT_TYPE mn, vector<MAT_TYPE> &vr, OUT MAT_TYPE & p);
void decimate(CDGMulDimData & dat, VECTOR_TYPE fwhm, CDGMulDimData &dout);
void fmgN(CDGMulDimData &A, CDGMulDimData &b, double p[8], double* p4, CDGMulDimData& rout);
void InvDef(CDGMulDimData &y, VECTOR_TYPE &d1, MAT_TYPE &M1_, MAT_TYPE &M2_, CDGMulDimData & invY);
void push(mwSize dm[], mwSize m, mwSize n, float def[], float pf[], float po[], float so[]);
MAT_TYPE BSplines(double *c, DimVec dim, MAT_TYPE &x0, MAT_TYPE &x1, MAT_TYPE &x2, DoubleVec &deg);
void defs(CDGMulDimData & TWrap_d, int z, MAT_TYPE &x0, MAT_TYPE &y0, VECTOR_TYPE &z0, MAT_TYPE &M, sMask & pmsk,
	OUT VECTOR_TYPE &x1, OUT VECTOR_TYPE &y1, OUT VECTOR_TYPE &z1);
MAT_TYPE transf(MAT_TYPE & B1, MAT_TYPE & B2, VECTOR_TYPE B3, Eigen::Map<MAT_TYPE> & T);
MAT_TYPE transf(MAT_TYPE & B1, MAT_TYPE & B2, VECTOR_TYPE B3, CDGMulDimData & T);
//Estimate cluster parameters			
void EstimateHistrogramParameters(IN OUT VECTOR_TYPE &WrapPrior, IN OUT ChanVec& chans, IN OUT CVariableLL &ll, IN sEstimateParam & eparam,
	IN MAT_TYPE &vr0, IN std::vector<sWrapBUF> &WrapBufs);
void EstimateClusterParameters(IN OUT VECTOR_TYPE &WrapPrior, IN OUT vector<MAT_TYPE> &vr, IN OUT VECTOR_TYPE & mg, IN OUT  MAT_TYPE & mn, IN OUT CVariableLL &ll,
	IN sEstimateParam & eparam, IN VECTOR_TYPE &lkp, IN std::vector<sWrapBUF> &WrapBufs, IN MAT_TYPE &vr0, int N, int &iter, int &iter1);
void InitHistrogram(IN OUT ChanVec& chans, IN sEstimateParam & eparam, IN std::vector<sWrapBUF> &WrapBufs);
void InitGaussMix(IN OUT vector<MAT_TYPE> &vr, IN OUT VECTOR_TYPE & mg, IN OUT  MAT_TYPE & mn, IN sEstimateParam & eparam, IN VECTOR_TYPE &lkp,
	IN std::vector<sWrapBUF> &WrapBufs, IN MAT_TYPE &vr0, int N);
void LoadChans(IN sChannels& channels, IN sRegData regdata, OUT ChanVec &chans);
void LoadWrapBufs(IN sChannels& channels, IN sRegData regdata, IN MAT_TYPE &M, IN ChanVec &chans, IN sTPMPriors & tpm, IN int objsample, OUT std::vector<sWrapBUF> &WrapBufs,
	IN OUT sEstimateParam & eparam, OUT MAT_TYPE &vr0);
void InitBfll(IN ChanVec &chans, IN OUT std::vector<sWrapBUF> &WrapBufs, OUT double &llrb);
void ObjectiveAndDerivateMog(OUT MAT_TYPE &wt1, OUT  MAT_TYPE &wt2, IN OUT vector<MAT_TYPE> &vr, IN OUT VECTOR_TYPE &WrapPrior,
	IN OUT VECTOR_TYPE & mg, IN OUT  MAT_TYPE & mn, IN sEstimateParam & eparam, IN VECTOR_TYPE &lkp, IN sWrapBUF&buf,IN vector<MAT_TYPE> &pr, IN DimVec &d, IN int n);
void ObjectiveAndDerivateNoMog(OUT MAT_TYPE &wt1, OUT  MAT_TYPE &wt2, IN OUT VECTOR_TYPE &WrapPrior,
	IN ChanVec&chans, IN sEstimateParam & eparam, IN DimVec &d, IN sWrapBUF&buf, IN int n);
void LineSearch(IN ChanVec & chans, int n, IN MAT_TYPE & Alpha, IN MAT_TYPE & Beta, vector<sWrapBUF> & WrapBufs,
	IN sEstimateParam & eparam, IN VECTOR_TYPE &lkp, bool busemog, IN OUT VECTOR_TYPE &WrapPrior,
	IN OUT vector<MAT_TYPE> &vr, IN OUT VECTOR_TYPE & mg, IN OUT  MAT_TYPE & mn, IN OUT CVariableLL &ll);
VECTOR_TYPE RandnVec(int N);
void UpdateAlphaBeta(MAT_TYPE & Alpha, MAT_TYPE & Beta, sChan & chan, MAT_TYPE &wt1, MAT_TYPE &wt2, int z);
void ReFillMix(IN OUT vector<MAT_TYPE> &vr, IN OUT VECTOR_TYPE & mg, IN OUT  MAT_TYPE & mn, IN VECTOR_TYPE &lkp, IN sEstimateParam & eparam, int N);
void EstimateDeformationsMog(IN sEstimateParam & eparam, IN OUT VECTOR_TYPE & mg, IN OUT  MAT_TYPE & mn, vector<MAT_TYPE> & vr,
	IN std::vector<sWrapBUF> &WrapBufs, IN OUT VECTOR_TYPE &WrapPrior, IN VECTOR_TYPE &lkp, IN OUT CVariableLL & ll);
void EstimateDeformationsNoMog(IN sEstimateParam & eparam, ChanVec& chans,
	IN std::vector<sWrapBUF> &WrapBufs, IN OUT VECTOR_TYPE &WrapPrior, IN OUT CVariableLL & ll);
void UpdateTWrap(IN OUT CVariableLL & ll, IN MAT_TYPE & x0, IN MAT_TYPE & y0, IN VECTOR_TYPE & z0, IN vector<sWrapBUF> & WrapBufs,
	IN VECTOR_TYPE & param, IN OUT sEstimateParam & eparam, IN VECTOR_TYPE &WrapPrior, IN CDGMulDimData & TWrap,
	IN sTPMPriors & tpm, IN MAT_TYPE MT, IN MAT_TYPE M, int iter, bool bDecreaseParamByIter, IN VECTOR_TYPE &sk);
float ComputeObjFunction(VECTOR_TYPE & sk, CDGMulDimData & TWrap, IN VECTOR_TYPE param);
void Fmg(CDGMulDimData & alpha, CDGMulDimData & beta, double _param[10], CDGMulDimData &X);
void Vel2Mom(double _param[8], CDGMulDimData & datain, CDGMulDimData & dataout);
float ComputeObjFunction(VECTOR_TYPE & sk, CDGMulDimData & TWrap, IN VECTOR_TYPE param);
void defs(CDGMulDimData & TWrap_f, int z, MAT_TYPE &x0, MAT_TYPE &y0, VECTOR_TYPE &z0, MAT_TYPE &M, sMask & pmsk,
	OUT VECTOR_TYPE &x1, OUT VECTOR_TYPE &y1, OUT VECTOR_TYPE &z1);
void latentnopar(IN vector<VECTOR_TYPE_F> &f, IN vector<VECTOR_TYPE_F> &bf, ChanVec& chans, IN  MAT_TYPE_F B, IN  VECTOR_TYPE& wp, OUT MAT_TYPE_F &Q, OUT double &ll);
void latent(IN  vector<VECTOR_TYPE_F> &f, IN  vector<VECTOR_TYPE_F> &bf, IN  VECTOR_TYPE & mg, IN  MAT_TYPE& mn, IN  vector<MAT_TYPE>& vr,
	IN VECTOR_TYPE &lkp, IN  MAT_TYPE B, IN  VECTOR_TYPE& wp, OUT MAT_TYPE &Q, OUT double &ll);
MAT_TYPE transf(MAT_TYPE & B1, MAT_TYPE & B2, VECTOR_TYPE B3, Eigen::Map<MAT_TYPE> & T);
MAT_TYPE transf(MAT_TYPE & B1, MAT_TYPE & B2, VECTOR_TYPE B3, CDGMulDimData & T);
void LoadBuffer(sTPMPriors & tpmpriors, sitk::Image & MovImg, sRegData & regdata, int nsample);
MAT_TYPE AffReg(sRegData & regdata, sTPMPriors & priors, MAT_TYPE M);
MAT_TYPE AffinePass(sitk::Image & objimage, int objsample, double fwhm, sTPMPriors & priors, MAT_TYPE affinemat, PriorsType priorstype);
VECTOR_TYPE VoxelSize(MAT_TYPE & mat);
void AffinePriors(PriorsType type, OUT VECTOR_TYPE &mu, OUT MAT_TYPE & isigma);
VECTOR_TYPE Matrix2Param(MAT_TYPE &M);
MAT_TYPE Param2Matrix(VECTOR_TYPE & param);
sitk::AffineTransform MatToAffine(MAT_TYPE & mat);
VECTOR_TYPE SmoothKernel(double fwhm, int xstart, int xend, KernelSmoothMode mode, double eps);
MAT_TYPE GetMatFromImageDirOriSpace(DoubleVec Direction, DoubleVec Origin, DoubleVec & space, bool bApplyMatlabIndexOffset, bool bITKSpaceToNifitySpace);
void ITKSpaceToNifitySpace(DoubleVec & Direction, DoubleVec & Origin);
void ScaleFactor(double & sf0, double & sf1, float pmax, float pmin, double cmin, double cmax);
VECTOR_TYPE SkFromVx(VECTOR_TYPE & vx, int nsample);
int BSplinec(DoubleVec &degs, double* volin, DimVec & dim, double * c);
VECTOR_TYPE BSplines(double *c, DimVec dim, VECTOR_TYPE &x0, VECTOR_TYPE &x1, VECTOR_TYPE &x2, DoubleVec &deg);
int BSplinesD(double *c, DimVec dim, DoubleVec &deg, VECTOR_TYPE &x1, VECTOR_TYPE &x2, VECTOR_TYPE &x3,
	OUT VECTOR_TYPE&f, VECTOR_TYPE &d1s, VECTOR_TYPE &d2s, VECTOR_TYPE &d3s);
void ApplyMatlabIndexOffset(MAT_TYPE & mat);
bool WrapPass(sTpmInputs& TInput, sTPMPriors & tpm, MAT_TYPE affinemat, sWrapResults &wresult);
void TpmWrite8(sTpmInputs & TInput, sTPMPriors & tpm, IN sWrapResults & res, sWriteImgs & sWrites);
MAT_TYPE Mat4x4FromAffineParam(DoubleVec & dp);
MAT_TYPE EigenMatFromDoubleVec(DoubleVec & dp, int row, int column);
MAT_TYPE_F AsSingle(MAT_TYPE & m);
VECTOR_TYPE_F AsSingle(VECTOR_TYPE & v);
MAT_TYPE AsDouble(const MAT_TYPE_F & m);
VECTOR_TYPE AsDouble(const VECTOR_TYPE_F & v);
MAT_TYPE Sqrtm(MAT_TYPE m);
VECTOR_TYPE Sqrtm(VECTOR_TYPE m);
MAT_TYPE Expm(MAT_TYPE m);
MAT_TYPE Logm(MAT_TYPE m);