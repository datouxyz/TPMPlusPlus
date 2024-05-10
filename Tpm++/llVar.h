#pragma once

class CVariableLL
{
public:
	CVariableLL()
	{

	}
	double GetLL()
	{
		return m_bDouble ? m_lldouble : m_llsingle;
	}
	void DoubleToSingle()
	{
		assert(m_bDouble);
		m_bDouble = false;
		m_llsingle = m_lldouble;
	}
	void SingleToDouble()
	{
		assert(!m_bDouble);
		m_bDouble = true;
		m_lldouble = m_llsingle;
	}
	bool bDouble() { return m_bDouble; }
	void IncreaseLL(double inc)
	{
		if (m_bDouble)
		{
			m_lldouble += inc;
		}
		else
		{
			m_llsingle += inc;
		}
	}
	void SetLL(double newll)
	{
		if (m_bDouble)
			m_lldouble = newll;
		else m_llsingle = newll;
	}
protected:
	bool m_bDouble{true};
	double m_lldouble{0.0};
	float m_llsingle{ 0.0f};
};