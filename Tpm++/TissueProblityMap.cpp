#include "stdafx.h"

#include "ParamaterDef.h"
#include "TissueProblityMap.h"


#include "EigenMatFunc.h"
#include "TPMInternalFuncs.h"
VECTOR_TYPE VoxelSize(MAT_TYPE & mat)
{	
	
	return mat.block(0, 0, 3, 3).array().pow(2.0).rowwise().sum().sqrt();
	//MAT_TYPE m33 = mat.block(0,0,3,3);

	//MAT_TYPE p2 = m33.array()*m33.array();
	//VECTOR_TYPE ret(3); ret << sqrt(p2.row(0).sum()), sqrt(p2.row(1).sum()), sqrt(p2.row(2).sum());
	/*cout << "ret" << ret;
	ret[0] = sqrt(mat(0,0) * mat(0,0)) + sqrt(mat(0,1) * mat(0,1)) + sqrt(mat(0,2) * mat(0,2));
	ret[1] = sqrt(mat(1,0) * mat(1,0)) + sqrt(mat(1,1) * mat(1,1)) + sqrt(mat(1,2) * mat(1,2));
	ret[2] = sqrt(mat(2,0) * mat(2,0)) + sqrt(mat(2,1) * mat(2,1)) + sqrt(mat(2,2) * mat(2,2));
	cout << "ret" << ret;*/
	//return ret;
}


MAT_TYPE GetMatFromImageDirOriSpace(DoubleVec  dir, DoubleVec  ori, DoubleVec & space, bool bApplyMatlabIndexOffset,bool bITKSpaceToNifitySpace)
{
	/*DoubleVec dir = image.GetDirection();
	DoubleVec ori = image.GetOrigin();
	DoubleVec space = image.GetSpacing();
*/
	if (bITKSpaceToNifitySpace)
	{
		ITKSpaceToNifitySpace(dir, ori);
	}
	

	VECTOR_TYPE vori = ToVecType(ori);
	MAT_TYPE MATRot(3, 3);
	MATRot << dir[0], dir[1], dir[2],
		dir[3], dir[4], dir[5],
		dir[6], dir[7], dir[8];

	MAT_TYPE MATSpace(3, 3); MATSpace << space[0], 0, 0,
		0, space[1], 0,
		0, 0, space[2];

	MAT_TYPE ret(4, 4);
	ret <<
		dir[0] * space[0], dir[1] * space[1], dir[2] * space[2], ori[0],
		dir[3] * space[0], dir[4] * space[1], dir[5] * space[2], ori[1],
		dir[6] * space[0], dir[7] * space[1], dir[8] * space[2], ori[2],
		0, 0, 0, 1;

	if (bApplyMatlabIndexOffset)
	{
		ApplyMatlabIndexOffset(ret);
	}	
	//cout << endl<< ret<<endl;
	return ret;
}

void ITKSpaceToNifitySpace(DoubleVec & Direction, DoubleVec & Origin)
{
	for (int i = 0; i < 6; i++)
	{
		Direction[i] *= -1;
	}

	{
		Origin[0] *= -1;
		Origin[1] *= -1;
	}
}

sTpmVol3D::sTpmVol3D(sitk::Image & tpmimg4d)
{
	m_Ori.resize(3);
	m_Dir.resize(9);
	m_Spacing.resize(3);

	m_Ori4D = tpmimg4d.GetOrigin();
	m_Ori[0] = m_Ori4D[0];
	m_Ori[1] = m_Ori4D[1];
	m_Ori[2] = m_Ori4D[2];

	m_Spacing4D = tpmimg4d.GetSpacing();
	m_Spacing[0] = m_Spacing4D[0];
	m_Spacing[1] = m_Spacing4D[1];
	m_Spacing[2] = m_Spacing4D[2];

	m_Dir4D = tpmimg4d.GetDirection();
	m_Dir[0] = m_Dir4D[0];
	m_Dir[1] = m_Dir4D[1];
	m_Dir[2] = m_Dir4D[2];

	m_Dir[3] = m_Dir4D[4];
	m_Dir[4] = m_Dir4D[5];
	m_Dir[5] = m_Dir4D[6];

	m_Dir[6] = m_Dir4D[8];
	m_Dir[7] = m_Dir4D[9];
	m_Dir[8] = m_Dir4D[10];
}
void sTpmVol3D::SetToImage(sitk::Image & img)
{
	if (img.GetDimension() == 3)
	{
		img.SetDirection(m_Dir);
		img.SetOrigin(m_Ori);
		img.SetSpacing(m_Spacing);
	}
	else if(img.GetDimension() == 4)
	{
		img.SetDirection(m_Dir4D);
		img.SetOrigin(m_Ori4D);
		img.SetSpacing(m_Spacing4D);
	}	
}

void DirOriSpaceUpDim(IN DoubleVec  ODirs, IN DoubleVec  OOris, IN DoubleVec  OSpacings,
	OUT DoubleVec & Dirs, DoubleVec & Oris, DoubleVec & Spacings)
{
	
	Oris = OOris;
	Oris.push_back(0.0);

	Spacings = OSpacings;
	Spacings.push_back(1.0);

	int olddim = (int)OOris.size();
	int newdim = (int)OOris.size()+1;
	Dirs.resize(newdim*newdim);
	for(int irow = 0; irow< newdim; irow++)
	for (int icolumn = 0; icolumn < newdim; icolumn++)
	{
		double & d = Dirs[irow*newdim+icolumn];
		if (irow < olddim && icolumn<olddim)
		{
			d = ODirs[irow*olddim + icolumn];
		}
		else if(icolumn==irow)
		{
			d = 1.0;
		}
		else
		{
			d = 0.0;
		}
	}
}
void SplitDirOriSpacingOn3DVol(OUT DoubleVec & Dirs, DoubleVec & Oris, DoubleVec & Spacings,
	IN DoubleVec  ODirs, IN DoubleVec  OOris, IN DoubleVec  OSpacings)
{
	int nOriDim = (int)OOris.size();
	for (int idir = 0; idir < nOriDim; idir++)
		for (int iaxis = 0; iaxis < nOriDim; iaxis++)
		{
			if (idir <3 && iaxis <3) //(i%oridim.size() != m_iSplitDim)
			{
				Dirs.push_back(ODirs[idir*nOriDim + iaxis]);
			}
		}

	for (int i = 0; i < nOriDim; i++)
	{
		if (i <3)
		{
			Oris.push_back(OOris[i]);
			Spacings.push_back(OSpacings[i]);
		}
	}
}

void SplitDirOriSpacing(OUT DoubleVec & Dirs, DoubleVec & Oris, DoubleVec & Spacings,
	IN DoubleVec  ODirs, IN DoubleVec  OOris, IN DoubleVec  OSpacings, int iSplitDim)
{
	int nOriDim = (int)OOris.size();
	for (int idir = 0; idir < nOriDim; idir++)
		for (int iaxis = 0; iaxis < nOriDim; iaxis++)
		{
			if (idir != iSplitDim && iaxis != iSplitDim) //(i%oridim.size() != m_iSplitDim)
			{
				Dirs.push_back(ODirs[idir*nOriDim + iaxis]);
			}
		}

	for (int i = 0; i < nOriDim; i++)
	{
		if (i != iSplitDim)
		{
			Oris.push_back(OOris[i]);
			Spacings.push_back(OSpacings[i]);
		}
	}
}
MAT_TYPE sTpmVol3D::GetMat(bool bApplyMatlabIndexOffset,bool bItkSpaceToNifity)
{
	DoubleVec dir = m_Dir, ori = m_Ori;
	if (bItkSpaceToNifity)
	{
		ITKSpaceToNifitySpace(dir, ori);
	}
	return GetMatFromImageDirOriSpace(m_Dir, m_Ori, m_Spacing, bApplyMatlabIndexOffset, bItkSpaceToNifity);
}


bool SaveStringToFile(string s, string filepath);
string LoadStringFromFile(string filepath);
void LoadChans(sChannels & channels, IN sRegData regdata, OUT ChanVec &chans);
void TpmWrite8(sTpmInputs & TInput, sTPMPriors & tpm, IN sWrapResults & res, sWriteImgs & sWrites);

MAT_TYPE LDivide(MAT_TYPE &A, MAT_TYPE &B);
MAT_TYPE_F LDivide(MAT_TYPE_F &A, MAT_TYPE_F &B);

int GetNumPhyCpu();
int TestTPM()
{
	Eigen::initParallel();
	int ncpu = GetNumPhyCpu();

	Eigen::setNbThreads(ncpu);
	omp_set_num_threads(ncpu);

	sitk::Image Img0 = sitk::ReadImage("CSDN.nii");	
	sitk::Image TPM = sitk::ReadImage("E:\\TPM.nii", sitk::sitkFloat64);
	sTissues tss = {sTissue(1,true,false,true,true),  sTissue(1,true,true,false,false),
					sTissue(2,true,false,false,false),sTissue(3,true,false,false,false),
					sTissue(4,true,false,false,false),sTissue(2,false,false,false,false) };
	
	sChannels Channels(1);
	Channels[0].m_Img = Img0;
	Channels[0].m_biascorrect = true;
	Channels[0].m_biasfield = true;
	sTpmInputs tinput(tss, Channels);
	tinput.m_dFieldInverse= true ;
	tinput.m_dFieldForward= true ;
	tinput.m_Mrf= false ;
	tinput.m_Cleanup= true ;

	// std::string strJson;
	// bool bEncode = struct2x::JSONEncoder(strJson) << tinput;
	// assert(bEncode);
	// std::cout << "json is :" << strJson.c_str() << std::endl;
		
	// //std::vector<sTissue> tisses;
	
	// bool bDecode = struct2x::JSONDecoder(strJson.c_str()) >> tinput;
	// assert(bDecode);
	
	sTPMPriors tpm;
	tpm.LoadPriors(TPM);		
	{
		{
			sWrapResults r;		
			sRegData regdata(Channels[0].m_Img, tinput.m_nsample, tinput.m_fwhm);
			sWriteImgs sWrites;
			TpmWrite8(tinput, tpm, r, sWrites);
		}		
		MAT_TYPE Affine = tinput.m_InitAffine;
		{			
			Affine = AffinePass(tinput.m_Channels[0].m_Img, tinput.m_nsample,(tinput.m_fwhm+1)*16,tpm,Affine, tinput.m_ptype);
			Affine = AffinePass(tinput.m_Channels[0].m_Img, tinput.m_nsample, tinput.m_fwhm,	  tpm,Affine, tinput.m_ptype);
		}			
		sWrapResults wpres;
		WrapPass(tinput, tpm,Affine, wpres);
		// {
		// 	JSONNode root(JSON_NODE);
		// 	root.set_name("sWrapResults");
		// 	wpres.SaveToJson(root);
		// 	string jsonstr = root.write_formatted();

		// 	SaveStringToFile(jsonstr, "e:\\result1.json");
		// }
	}
	//}
	return 0;
}



void ndgrid_vec3(ByteVec & outdatas, int XBegin, int XEnd, int XStep,
	int YBegin, int YEnd, int YStep, int ZBegin, int ZEnd, int ZStep)
{
	int nColumn = (XEnd - XBegin) / (XStep);
	int nRow = (YEnd - YBegin) / (YStep);
	int nDepth = (ZEnd - ZBegin) / (ZStep);

	DimVec dim = { nColumn,nRow,nDepth};
	int nElenment = DimElCount(dim);
	outdatas.resize(sizeof(VERTEX3D)*nElenment);
	//DGDataInitilize(DGDataType_Vec3, outdatas, false, nElenment);
	
	VERTEX3D * pDataStart = (VERTEX3D*)&outdatas[0];

	for (int x =0;x<nColumn;x++)
	{
		for (int y = 0;y < nRow;y++)
		{
			for (int z = 0;z < nDepth;z++)
			{
				UINT iel = Dim3index(x,y,z,dim);
				VERTEX3D & pv = pDataStart[iel];
				pv.x() = (float)XBegin + XStep * x;
				pv.y() =(float) YBegin + YStep * y;
				pv.z() = (float)ZBegin + ZStep * z;
			}
		}
	}	
}
sTPMPriors::sTPMPriors()
:m_TPMMat(4,4)
{
	m_TPMMat = MAT_TYPE::Identity(4, 4);
}

int Slice0(double mat[], double image[], int xdim1, int ydim1, double vol[], int xdim2,int ydim2,int zdim2,double  background,double scale[],double offset[])
{
	double y, x2, y2, z2, s2, dx3 = mat[0], dy3 = mat[1], dz3 = mat[2], ds3 = mat[3];
	int t = 0;

	x2 = mat[12] + 0 * mat[8];
	y2 = mat[13] + 0 * mat[9];
	z2 = mat[14] + 0 * mat[10];
	s2 = mat[15] + 0 * mat[11];

	for (y = 0; y < ydim1; y++)
	{
		double x;
		double x3 = x2 + y * mat[4];
		double y3 = y2 + y * mat[5];
		double z3 = z2 + y * mat[6];
		double s3 = s2 + y * mat[7];

		for (x = 0; x < xdim1; x++)
		{
			int ix4, iy4, iz4;
			s3 += ds3;
			if (s3 == 0.0) return -1;
			
			ix4 = (int)floor(((x3 += dx3) / s3) -0.5);
			iy4 = (int)floor(((y3 += dy3) / s3));// -0.5);
			iz4 = (int)floor(((z3 += dz3) / s3));// -0.5);
			if (iz4 >= 0 && iz4 < zdim2 && iy4 >= 0 && iy4 < ydim2 && ix4 >= 0 && ix4 < xdim2)
			{
				image[t] = scale[iz4] * (double)(vol[iz4*xdim2*ydim2 + xdim2 * iy4+ ix4 ]) + offset[iz4];
			}
			else 
				image[t] = background;
			
			t++;
		}
	}
	return 0;
}


VECTOR_TYPE SkFromVx(VECTOR_TYPE & vx, int nsample)
{
	VECTOR_TYPE one3 = VECTOR_TYPE::Ones(3);
	VECTOR_TYPE SampleInVx = nsample * 1.0 / vx.array(); 
	VECTOR_TYPE rounds = vMax(one3, SampleInVx);
	Rounds(rounds);
	
	VECTOR_TYPE sk = vMax(one3, rounds);
	return sk;
}

double FugeFactor(double fwhm, VECTOR_TYPE &vx, VECTOR_TYPE &sk)
{
	//fudge factor
	//% Fudge Factor - to(approximately) account for
	//	% non - independence of voxels
	//	s = (fwhm + mean(vx)) / sqrt(8 * log(2));                 % Standard deviation
	//	ff = prod(4 * pi*(s. / vx. / sk). ^ 2 + 1) ^ (1 / 2);
	//Standard deviation	
	double s = (fwhm + vx.mean()) / sqrt(8 * log(2.0));

#define FFDIM(Dim) 4*PI*pow(s/vx[Dim]/sk[Dim],2)+1 

	DoubleVec fv3 = { FFDIM(0),FFDIM(1),FFDIM(2) };

	return sqrt(fv3[0] * fv3[1] * fv3[2]);
}


void sRegData::InitCoords(SizeVec movsize)
{
	int nx1 = (int)ceil(movsize[0] / m_sk[0]);
	int nx2 = (int)ceil(movsize[1] / m_sk[1]);
	int nx3 = (int)ceil(movsize[2] / m_sk[2]);

	x1 = MAT_TYPE::Zero(nx1, nx2);
	x2 = MAT_TYPE::Zero(nx1, nx2);
	x3 = VECTOR_TYPE::Zero(nx3);

	for (int ix1 = 0; ix1 < nx1; ix1++)
	{
		for (int ix2 = 0; ix2 < nx2; ix2++)
		{
			x1(ix1, ix2) = ix1 * m_sk[0] + 1;
			x2(ix1, ix2) = ix2 * m_sk[1] + 1;	
		}
	}	
	for (int ix3 = 0; ix3 < nx3; ix3++)
		x3(ix3) = ix3 * m_sk[2] + 1;
}
sitk::Image ReSampleImgAsFloat(sitk::Image & img,int nsample)
{
	SizeVec movsize = img.GetSize();
	
	MAT_TYPE MG = GetMatFromImageDirOriSpace(img.GetDirection(), img.GetOrigin(), img.GetSpacing(), true, true);//GetMatFromImageDirOriSpace(MovImg,true) ;
	
	VECTOR_TYPE vx = VoxelSize(MG);
	VECTOR_TYPE sk = SkFromVx(vx, nsample);

	unsigned int nx1 = (int)ceil(movsize[0] / sk[0]);
	unsigned int nx2 = (int)ceil(movsize[1] / sk[1]);
	unsigned int nx3 = (int)ceil(movsize[2] / sk[2]);
	
	DoubleVec Ori = img.GetOrigin();
	DoubleVec Space = img.GetSpacing();
	DoubleVec Dir = img.GetDirection();

	Space[0] *= sk[0];
	Space[1] *= sk[1];
	Space[2] *= sk[2];

	 SizeVec fsize ={nx1,nx2,nx3} ;
	 sitk::Image fimg = sitk::Resample(img, fsize, sitk::Transform(), sitk::sitkNearestNeighbor,
	 	Ori, Space, Dir, 0.0, sitk::sitkFloat32);
	return fimg;
}
sRegData::sRegData(sitk::Image & movimg, int nsample,double fwhm)
	: m_MG(4, 4)
{
	m_MG = GetMatFromImageDirOriSpace(movimg.GetDirection(), movimg.GetOrigin(), movimg.GetSpacing(), true, true);//GetMatFromImageDirOriSpace(MovImg,true) ;
	
	m_vx = VoxelSize(m_MG);
	m_sk = SkFromVx(m_vx, nsample);
	SizeVec movsize = movimg.GetSize();

	m_FudgeFactor = FugeFactor(fwhm, m_vx, m_sk);
	InitCoords(movsize);
}

void LoadBuffer(sTPMPriors & tpmpriors , sitk::Image & MovImg, sRegData & regdata, int nsample)
{		
	SizeVec movsize = MovImg.GetSize();		
	sitk::Image fimg = ReSampleImgAsFloat(MovImg,nsample);
	
	float* pfdata = fimg.GetBufferAsFloat();
	int ncomp = fimg.GetNumberOfComponentsPerPixel();//*fimg.GetSizeOfPixelComponent();
	double pmin = DBL_MAX, pmax = -DBL_MAX;
	for (int ipixel = 0; ipixel < fimg.GetNumberOfPixels(); ipixel++)
	{
		float & pixel = pfdata[ipixel*ncomp];
		pmin = min<double>(pixel, pmin);
		pmax = max<double>(pixel, pmax);
	}

	double sf0 = 0, sf1 = 0;

	ScaleFactor(sf0, sf1, (float)pmax, (float)pmin, 1, 4000);

	VECTOR_TYPE Historgram = VECTOR_TYPE::Zero(4000);

	SizeVec d = fimg.GetSize();
	UINT framesize = d[0] * d[1] * ncomp;

	int np = (int)fimg.GetNumberOfPixels();
	for (int ipixel = 0; ipixel < np; ipixel++)
	{
		float & pixel = pfdata[ipixel*ncomp];
		if (isnan(pixel) || isinf(pixel) || pixel == 0)
			continue;
		int  p =(int)round(pixel * sf0 + sf1)-1;
		Historgram[p]++;
	}

	VECTOR_TYPE hcums = CumSum(Historgram);
	double sum = Historgram.sum();

	for (auto &ch : hcums)
	{
		ch /= sum;
	}

	pmin = FLT_MAX, pmax = -FLT_MAX;
	double n_mn = 0, n_mx = 0;
	for (int i = 0; i < hcums.size(); i++)
	{
		double & cum = hcums[i];
		if (cum > 0.0005)
		{
			n_mn = i+1;
			break;
		}
	}
	for (int i = 0; i < hcums.size(); i++)
	{
		double & cum = hcums[i];
		if (cum > 0.9995)
		{
			n_mx = i + 1;
			break;
		}
	}

	n_mn = (n_mn - sf1) / sf0;
	n_mx = (n_mx - sf1) / sf0;

	double csf0 = 0;
	double csf1 = 0;
	ScaleFactor(csf0, csf1, (float)n_mx,(float) n_mn, 0, 255);
	
	regdata.m_Buffers.resize(d[2]);

	for (int iz = 0; iz < (int)d[2]; iz++)
	{
		float  *pframe = &pfdata[framesize*iz];
		sRegBuffer & srb = regdata.m_Buffers[iz];
		srb.m_g = VECTOR_TYPE(framesize);

		sMask &msk = srb.m_msk;
		msk = sMask(framesize);
		
		for (int ipixel = 0; ipixel < (int)framesize; ipixel++)
		{
			float  pixel = pframe[ipixel*ncomp];

			bool bSkip = isnan(pixel) || isinf(pixel) || pixel == 0;
			msk.SetMaskAndNM(ipixel,!bSkip);
		
			if (!bSkip)
			{			
				srb.m_g[ipixel] = clamp<double>(round(pixel * csf0 + csf1), 0, 255);
			}			
		}
		srb.m_g = srb.m_msk.MaskVector(srb.m_g); 		
	}
}

MAT_TYPE AffinePass(sitk::Image & objimage, int objsample, double fwhm, sTPMPriors & priors, MAT_TYPE affinemat, PriorsType priorstype)
{	
	sRegData regdata(objimage,objsample, fwhm);
	LoadBuffer(priors, objimage, regdata, objsample);
	MAT_TYPE M = AffReg( regdata, priors, affinemat);
	return M;
}


void ScaleFactor(double & sf0,double & sf1, float pmax , float pmin, double cmin, double cmax)
{
	MAT_TYPE m1(2, 2); m1 << pmin, 1, pmax, 1;
	VECTOR_TYPE r1(2); r1 << cmin, cmax;
	VECTOR_TYPE r = m1.inverse()*r1;
	sf0 = r(0);
	sf1 = r(1);
}
