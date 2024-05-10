#pragma once
enum PriorsType
{
	PriorsType_mni = 0,
	PriorsType_imni,
	PriorsType_rigid,
	PriorsType_subject,
	PriorsType_eastern,
	PriorsType_none = 0xffffffff
};


struct sTissue
{
public:
	sTissue(int ngaus, bool NativeSpace, bool DartelImported, bool Modulated, bool UnModulated)
		:m_ngaus(ngaus), m_NativeTissue_NativeSpace(NativeSpace), m_NativeTissue_DartelImported(DartelImported),
		m_WarpTissue_Modulated(Modulated), m_WarpTissue_UnModulated(UnModulated)
	{
	}
	sTissue()
	{

	}
	int m_ngaus{ 1 };
	//1
	bool m_NativeTissue_NativeSpace{ false };
	//2
	bool m_NativeTissue_DartelImported{ false };
	//3
	bool m_WarpTissue_Modulated{ false };
	//4
	bool m_WarpTissue_UnModulated{ false };
	template<typename T>
	void serialize(T& t) {
		SERIALIZE(t, m_ngaus, m_NativeTissue_NativeSpace, m_NativeTissue_DartelImported, m_WarpTissue_Modulated, m_WarpTissue_UnModulated);
	}
};
typedef std::vector<sTissue> sTissues;
struct sChannel
{
	sitk::Image m_Img;

	double m_Bisareg{ 1.0e-03 };
	double m_Biasfwhm{ 60 };
	bool m_biascorrect{ false };
	bool m_biasfield{ false };
	template<typename T>
	void serialize(T& t) {
		SERIALIZE(t, m_Bisareg, m_Biasfwhm, m_biascorrect, m_biasfield);
	}
};
typedef std::vector<sChannel> sChannels;


struct sTpmInputs
{
	sTpmInputs()
	{
		m_InitAffine = MAT_TYPE::Identity(4, 4);
	}
	sTpmInputs(sTissues &tss, sChannels &Channels)
		:m_reg{ 0, 0.001, 0.5, 0.05, 0.2 }, m_tissues(tss),
		m_bb(MAT_TYPE::Constant(NAN, 2, 3)), m_Channels(Channels)
	{
		m_lkp.clear();
		for (int k = 0, igauss = 0; k < tss.size(); k++)
		{
			auto & tiss = tss[k];
			for (int ig = 0; ig < tiss.m_ngaus; ig++)
			{
				m_lkp.push_back(k + 1);
			}
		}
		m_InitAffine = MAT_TYPE::Identity(4, 4);
	}
	bool AnyTissNaTNatSpace()
	{
		for (auto & ti : m_tissues)
		{
			if (ti.m_NativeTissue_NativeSpace)
				return true;
		}
		return false;
	}
	bool AnyTissNaTDartelSpace()
	{
		for (auto & ti : m_tissues)
		{
			if (ti.m_NativeTissue_DartelImported)
				return true;
		}
		return false;
	}

	bool AnyTissWarpModulated()
	{
		for (auto & ti : m_tissues)
		{
			if (ti.m_WarpTissue_Modulated)
				return true;
		}
		return false;
	}

	bool AnyTissWarpUnModulated()
	{
		for (auto & ti : m_tissues)
		{
			if (ti.m_WarpTissue_UnModulated)
				return true;
		}
		return false;
	}

	bool AnyChBiasCorrected()
	{
		for (auto & ch : m_Channels)
		{
			if (ch.m_biascorrect)
				return true;
		}
		return false;
	}
	bool AnyChBiasField()
	{
		for (auto & ch : m_Channels)
		{
			if (ch.m_biasfield)
				return true;
		}
		return false;
	}

	sChannels m_Channels;
	sTissues m_tissues;
	PriorsType m_ptype{ PriorsType_mni };

	DoubleVec m_reg;

	DoubleVec m_lkp;
	int m_nsample{ 3 };
	double m_fwhm{ 0 };
	MAT_TYPE m_bb;


	MAT_TYPE m_InitAffine;

	bool m_dFieldInverse{ true };
	bool m_dFieldForward{ true };
	bool m_Mrf{ true };
	int m_Cleanup{ 1 };//0,1,2

	VECTOR_TYPE m_WrapPrior;
	
	VECTOR_TYPE mg;
	MAT_TYPE mn;
	vector<MAT_TYPE> vr;
	CDGMulDimData m_TWrap;

	template<typename T>
	void serialize(T& t) {
		SERIALIZE(t, m_tissues, m_Channels, m_nsample, m_fwhm, m_reg, m_lkp, m_ptype,
			m_dFieldInverse, m_dFieldForward, m_Mrf, m_Cleanup);
	}

};