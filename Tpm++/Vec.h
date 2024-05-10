#pragma once
#include <vector>

template <typename T>
class Vec : public std::vector<T> {
public:
	Vec(unsigned __int64 size)
		:std::vector<T>(size) {}
	Vec(std::vector<T> & cp)
		: std::vector<T>(cp) {}

	T GetMax() const {
		if (empty()) throw std::invalid_argument("Vector is empty");
		return *std::max_element(begin(), end());
	}

	Vec operator+(const T& value) const {
		Vec result(size());
		for (size_t i = 0; i < size(); i++) {
			result[i] = at(i) + value;
		}
		return result;
	}
	template <typename U>
	Vec<U> castTo() const {
		Vec<U> result(size());
		for (size_t i = 0; i < size(); i++) {
			result[i] = static_cast<U>(at(i));
		}
		return result;
	}
	Vec operator-(const T& value) const {
		Vec result(size());
		for (size_t i = 0; i < size(); i++) {
			result[i] = at(i) - value;
		}
		return result;
	}

	Vec operator*(const T& value) const {
		Vec result(size());
		for (size_t i = 0; i < size(); i++) {
			result[i] = at(i) * value;
		}
		return result;
	}

	Vec operator/(const T& value) const {
		if (value == 0) throw std::invalid_argument("Division by zero");
		Vec result(size());
		for (size_t i = 0; i < size(); i++) {
			result[i] = at(i) / value;
		}
		return result;
	}
	Vec operator+(const Vec& b) const {
		if (size() != b.size()) throw std::invalid_argument("Vectors must be the same size");
		Vec result(size());
		for (size_t i = 0; i < size(); i++) {
			result[i] = at(i) + b[i];
		}
		return result;
	}

	Vec operator-(const Vec& b) const {
		if (size() != b.size()) throw std::invalid_argument("Vectors must be the same size");
		Vec result(size());
		for (size_t i = 0; i < size(); i++) {
			result[i] = at(i) - b[i];
		}
		return result;
	}

	Vec operator*(const Vec& b) const {
		if (size() != b.size()) throw std::invalid_argument("Vectors must be the same size");
		Vec result(size());
		for (size_t i = 0; i < size(); i++) {
			result[i] = at(i) * b[i];
		}
		return result;
	}

	Vec operator/(const Vec& b) const {
		if (size() != b.size()) throw std::invalid_argument("Vectors must be the same size");
		Vec result(size());
		for (size_t i = 0; i < size(); i++) {
			if (b[i] == 0) throw std::invalid_argument("Division by zero");
			result[i] = at(i) / b[i];
		}
		return result;
	}
};