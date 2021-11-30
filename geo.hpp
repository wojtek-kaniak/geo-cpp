#ifndef GEO_HPP_INCLUDED
#define GEO_HPP_INCLUDED
#include <array>
#include <optional>
#include <utility>
#include <iostream>
#include <cassert>

namespace geo
{

template <typename T>
class Matrix3x3
{
public:
	Matrix3x3(std::initializer_list<T> values)
	{
		assert(values.size() == vals.size());
		std::copy(values.begin(), values.end(), vals.begin());
	}

	T& operator[](std::pair<int, int> ix) noexcept
	{
		return vals[ix.second * 3 + ix.first];
	}

	const T& operator[](std::pair<int, int> ix) const noexcept
	{
		return vals[ix.second * 3 + ix.first];
	}

private:
	std::array<T, 9> vals;
};

template <typename T>
class Matrix2x2
{
public:
	Matrix2x2(std::initializer_list<T> values)
	{
		assert(values.size() == vals.size());
		std::copy(values.begin(), values.end(), vals.begin());
	}

	T& operator[](std::pair<int, int> ix) noexcept
	{
		return vals[ix.second * 2 + ix.first];
	}

	const T& operator[](std::pair<int, int> ix) const noexcept
	{
		return vals[ix.second * 2 + ix.first];
	}

private:
	std::array<T, 4> vals;
};

template <typename T>
class Point
{
public:
	constexpr Point(const T& x, const T& y) noexcept : x(x), y(y) {}
	T x, y;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Point<T>& point)
{
	os << '(' << point.x << ';' << point.y << ')';
	return os;
}

template <typename T>
static T gcd(T a, T b)
{
	while (b)
	{
		T c = b;
		b = a % b;
		a = c;
	}
	return a;
}

template <typename T>
class Fraction
{
public:
	constexpr Fraction(const T& numerator, const T& denominator)
	{
		num = numerator;
		den = (denominator != (T)0) ? denominator : throw std::domain_error("division by 0");
	}

	constexpr explicit operator double() const
	{
		return (double)num / (double)den;
	}

	constexpr Fraction<T> reduce() const noexcept
	{
		T div = gcd(num, den);
		return {num / div, den / div};
	}
	
	constexpr void reduce() noexcept
	{
		T div = gcd(num, den);
		num /= div;
		den /= div;
	}

	T num, den;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Fraction<T>& frac)
{
	os << frac.num << '/' << frac.den;
	return os;
}

template <typename T>
using pos_t = Point<T>;
template <typename T>
using seg_t = std::pair<pos_t<T>, pos_t<T>>;
template <typename T>
using frac_t = Fraction<T>;

// Determinant of a 3x3 matrix
template <typename T>
T det(Matrix3x3<T>& matrix) noexcept
{
	return
		+ matrix[{0, 0}] * matrix[{1, 1}] * matrix[{2, 2}]
		+ matrix[{0, 1}] * matrix[{1, 2}] * matrix[{2, 0}]
		+ matrix[{0, 2}] * matrix[{1, 0}] * matrix[{2, 1}]
		- matrix[{2, 0}] * matrix[{1, 1}] * matrix[{0, 2}]
		- matrix[{2, 1}] * matrix[{1, 2}] * matrix[{0, 0}]
		- matrix[{2, 2}] * matrix[{1, 0}] * matrix[{0, 1}];
}

// Determinant of a 2x2 matrix
template <typename T>
constexpr T det(Matrix2x2<T>& matrix) noexcept
{
	return matrix[{0, 0}] * matrix[{1, 1}] - matrix[{1, 0}] * matrix[{0, 1}];
}

template <typename T>
int sgn(T value) noexcept
{
	return (T(0) < value) - (value < T(0));
}

// Determine on which side of line seg given point is
// -1 - left
//  0 - collinear
//  1 - right
template <typename T>
int side(seg_t<T> seg, pos_t<T> point) noexcept
{
	Matrix3x3<T> matrix = {
		seg.first.x, seg.first.y, (T)1,
		seg.second.x, seg.second.y, (T)1,
		point.x, point.y, (T)1
	};

	return sgn<T>(det<T>(matrix));
}

// Determine whether both points lie on the same side of line seg
template <typename T>
bool same_side(seg_t<T> seg, pos_t<T> point1, pos_t<T> point2) noexcept
{
	return side<T>(seg, point1) == side<T>(seg, point2);
}

// Determine whether segment seg cointains given point
template <typename T>
bool seg_contains(seg_t<T> seg, pos_t<T> point) noexcept
{
	return side<T>(seg, point) == (T)0
		&& std::min(seg.first.x, seg.second.x) <= point.x <= std::max(seg.first.x, seg.second.x)
		&& std::min(seg.first.y, seg.second.y) <= point.y <= std::max(seg.first.y, seg.second.y);
}

// Determine whether given segments intersect
template <typename T>
bool seg_intersects(seg_t<T> seg1, seg_t<T> seg2) noexcept
{
	return
		(!same_side(seg1, seg2.first, seg2.second) && !same_side(seg2, seg1.first, seg1.second))
		|| seg_contains(seg1, seg2.first) || seg_contains(seg1, seg2.second)
		|| seg_contains(seg2, seg1.first) || seg_contains(seg2, seg1.second);
}

// Calculate the exact intersection point of given lines
template <typename T>
std::optional<pos_t<frac_t<T>>> seg_intersection(seg_t<T> seg1, seg_t<T> seg2) noexcept
{
	Matrix2x2 den = {
		seg1.first.x - seg1.second.x, seg1.first.y - seg1.second.y,
		seg2.first.x - seg2.second.x, seg2.first.y - seg2.second.y,
	};

	// determinant of the denominator
	T den_det = det<T>(den);

	// parallel lines
	if (den_det == (T)0)
		return {};

	Matrix2x2 m1 = {
		seg1.first.x, seg1.first.y,
		seg1.second.x, seg1.second.y
	};

	Matrix2x2 m2 = {
		seg2.first.x, seg2.first.y,
		seg2.second.x, seg2.second.y
	};

	Matrix2x2 xnum = {
		det<T>(m1), seg1.first.x - seg1.second.x,
		det<T>(m2), seg2.first.x - seg2.second.x,
	};
	
	Matrix2x2 ynum = {
		det<T>(m1), seg1.first.y - seg1.second.y,
		det<T>(m2), seg2.first.y - seg2.second.y,
	};

	pos_t<frac_t<T>> frac = {
		{det<T>(xnum), den_det},
		{det<T>(ynum), den_det}
	};

	frac.x.reduce();
	frac.y.reduce();

	return frac;
}
}

#endif
