defmodule ComplexTest do
  use ExUnit.Case
  doctest Complex

  test "Creation" do
    assert_close Complex.new(1.0, 1.0), %Complex{re: 1.0, im: 1.0}
    assert_close Complex.new(1.0, 2.0), %Complex{re: 1.0, im: 2.0}
    assert_close Complex.new(2.0), %Complex{re: 2.0, im: 0.0}

    a = Complex.new(1.0, 2.0)
    {r, theta} = Complex.to_polar(a)

    assert_close r, 2.23606797749979
    assert_close theta, 1.1071487177940904

    assert {:infinity, 0} == Complex.to_polar(:infinity)
    assert {:infinity, :math.pi()} == Complex.to_polar(:neg_infinity)
    assert {:nan, :nan} == Complex.to_polar(:nan)

    assert_close Complex.from_polar(1.0, :math.pi() / 2), %Complex{
      re: 6.123233995736766e-17,
      im: 1.0
    }

    assert Complex.from_polar({:infinity, :neg_infinity}) == %Complex{re: :nan, im: :nan}
    assert Complex.from_polar({:infinity, 0}) == %Complex{re: :infinity, im: 0}
    assert Complex.from_polar({:infinity, :math.pi()}) == %Complex{re: :neg_infinity, im: 0}

    assert Complex.from_polar({:infinity, :math.pi() / 2}) == %Complex{re: 0, im: :infinity}
    assert Complex.from_polar({:infinity, -:math.pi() / 2}) == %Complex{re: 0, im: :neg_infinity}
  end

  test "Parse from string" do
    assert_close Complex.parse("1.0+1.0i") |> elem(0), %Complex{re: 1.0, im: 1.0}
    assert_close Complex.parse("-1+1.0i") |> elem(0), %Complex{re: -1.0, im: 1.0}
    assert_close Complex.parse("1.+1.0i") |> elem(0), %Complex{re: 1.0, im: 1.0}
    assert_close Complex.parse("1-1.0i") |> elem(0), %Complex{re: 1.0, im: -1.0}

    tail = "12345"
    assert {%Complex{re: :infinity, im: :nan}, tail} == Complex.parse("+Inf-NaNi" <> tail)

    assert {%Complex{re: :neg_infinity, im: :infinity}, tail} ==
             Complex.parse("-Inf+Infi" <> tail)

    assert {%Complex{re: :infinity, im: :neg_infinity}, tail} == Complex.parse("Inf-Infi" <> tail)
    assert {%Complex{re: :nan, im: :nan}, tail} == Complex.parse("NaN+NaNi" <> tail)

    assert {%Complex{re: :infinity, im: 0}, tail} == Complex.parse("Inf" <> tail)
    assert {%Complex{re: :neg_infinity, im: 0}, tail} == Complex.parse("-Inf" <> tail)
    assert {%Complex{re: :nan, im: 0}, tail} == Complex.parse("NaN" <> tail)

    assert {%Complex{re: 0, im: :infinity}, ""} == Complex.parse("Infi")
    assert {%Complex{re: 0, im: :infinity}, ""} == Complex.parse("+Infi")
    assert {%Complex{re: 0, im: :neg_infinity}, ""} == Complex.parse("-Infi")
    assert {%Complex{re: 0, im: :nan}, ""} == Complex.parse("NaNi")
    assert {%Complex{re: 0, im: :nan}, ""} == Complex.parse("+NaNi")
    assert {%Complex{re: 0, im: :nan}, ""} == Complex.parse("-NaNi")

    assert :error == Complex.parse("123")
    assert_close Complex.parse("1+1i") |> elem(0), %Complex{re: 1.0, im: 1.0}
  end

  test "Arithmetic" do
    a = Complex.new(1.0, 2.0)
    b = Complex.new(3.0, 4.0)
    assert_close Complex.add(a, b), %Complex{re: 4.0, im: 6.0}
    assert_close Complex.subtract(a, b), %Complex{re: -2.0, im: -2.0}
    assert_close Complex.multiply(a, b), %Complex{re: -5.0, im: 10.0}
    assert_close Complex.divide(a, b), %Complex{re: 0.44, im: 0.08}
    assert_close Complex.divide(b, a), %Complex{re: 2.2, im: -0.4}
  end

  test "Arithmetic (no upcasts)" do
    a = 2.0
    b = 3.0
    assert_close Complex.add(a, b), 5
    assert_close Complex.subtract(a, b), -1
    assert_close Complex.multiply(a, b), 6
    assert_close Complex.divide(a, b), 0.666666
    assert_close Complex.divide(b, a), 1.5
  end

  test "add and subtract (non finite)" do
    for x <- [:nan, 1, -1, 0, :infinity, :neg_infinity] do
      assert Complex.add(:nan, x) == :nan
      assert Complex.add(x, :nan) == :nan
      assert Complex.subtract(:nan, x) == :nan
      assert Complex.subtract(x, :nan) == :nan
    end

    assert Complex.add(:infinity, :infinity) == :infinity
    assert Complex.add(:infinity, :neg_infinity) == :nan
    assert Complex.add(:neg_infinity, :infinity) == :nan
    assert Complex.add(:neg_infinity, :neg_infinity) == :neg_infinity

    assert Complex.subtract(:infinity, :infinity) == :nan
    assert Complex.subtract(:infinity, :neg_infinity) == :infinity
    assert Complex.subtract(:neg_infinity, :infinity) == :neg_infinity
    assert Complex.subtract(:neg_infinity, :neg_infinity) == :nan

    for x <- [-1, 0, 1] do
      assert Complex.add(:infinity, x) == :infinity
      assert Complex.add(:neg_infinity, x) == :neg_infinity
      assert Complex.subtract(:infinity, x) == :infinity
      assert Complex.subtract(:neg_infinity, x) == :neg_infinity
    end

    # property test
    values = [:infinity, :neg_infinity, :nan, -1, 0, 1]

    for re1 <- values, im1 <- values, re2 <- values, im2 <- values do
      assert Complex.add(Complex.new(re1, im1), Complex.new(re2, im2)) ==
               Complex.new(Complex.add(re1, re2), Complex.add(im1, im2))

      assert Complex.subtract(Complex.new(re1, im1), Complex.new(re2, im2)) ==
               Complex.new(Complex.subtract(re1, re2), Complex.subtract(im1, im2))
    end
  end

  test "multiply (non finite)" do
    assert Complex.multiply(:infinity, -1) == :neg_infinity
    assert Complex.multiply(:infinity, 0) == :nan
    assert Complex.multiply(:infinity, 1) == :infinity
    assert Complex.multiply(-1, :infinity) == :neg_infinity
    assert Complex.multiply(0, :infinity) == :nan
    assert Complex.multiply(1, :infinity) == :infinity

    assert Complex.multiply(:neg_infinity, -1) == :infinity
    assert Complex.multiply(:neg_infinity, 0) == :nan
    assert Complex.multiply(:neg_infinity, 1) == :neg_infinity
    assert Complex.multiply(-1, :neg_infinity) == :infinity
    assert Complex.multiply(0, :neg_infinity) == :nan
    assert Complex.multiply(1, :neg_infinity) == :neg_infinity

    for x <- [:infinity, :neg_infinity, -1, 0, 1] do
      assert Complex.multiply(:nan, x) == :nan
      assert Complex.multiply(x, :nan) == :nan
    end

    for x <- [:infinity, :neg_infinity], y <- [:infinity, :neg_infinity] do
      if x == y do
        assert Complex.multiply(x, y) == :infinity
      else
        assert Complex.multiply(x, y) == :neg_infinity
      end
    end

    assert Complex.multiply(:infinity, Complex.new(1, 2)) == Complex.new(:infinity, :infinity)

    assert Complex.multiply(:infinity, Complex.new(-3, -4)) ==
             Complex.new(:neg_infinity, :neg_infinity)

    assert Complex.multiply(:neg_infinity, Complex.new(1, 2)) ==
             Complex.new(:neg_infinity, :neg_infinity)

    assert Complex.multiply(:neg_infinity, Complex.new(-3, -4)) ==
             Complex.new(:infinity, :infinity)

    assert Complex.multiply(Complex.new(:infinity), Complex.new(:neg_infinity, 1)) ==
             Complex.new(:neg_infinity, :nan)

    assert Complex.square(Complex.new(1, :infinity)) == Complex.new(:neg_infinity, :infinity)
    assert Complex.square(Complex.new(:infinity, 1)) == Complex.new(:infinity, :infinity)
  end

  test "divide (non finite)" do
    for x <- [:infinity, :neg_infinity, :nan], y <- [:infinity, :neg_infinity, :nan] do
      assert Complex.divide(x, y) == :nan
      assert Complex.divide(x, 0) == x
      assert Complex.divide(x, 1) == x

      if y == :nan do
        assert Complex.divide(Complex.new(x), Complex.new(y)) == Complex.new(:nan, :nan)
      else
        assert Complex.divide(Complex.new(x), Complex.new(y)) == Complex.new(:nan)
      end

      assert Complex.divide(Complex.new(x), 0) == Complex.new(x, :nan)
      assert Complex.divide(Complex.new(x), 1) == Complex.new(x)
    end

    assert Complex.divide(:nan, 1) == :nan
    assert Complex.divide(1, :nan) == :nan

    assert Complex.divide(Complex.new(:nan), 1) == Complex.new(:nan, 0)
    assert Complex.divide(1, Complex.new(:nan)) == Complex.new(:nan, :nan)

    assert Complex.divide(:infinity, -1) == :neg_infinity
    assert Complex.divide(:neg_infinity, -1) == :infinity

    assert Complex.divide(Complex.new(:infinity), -1) == Complex.new(:neg_infinity, 0)
    assert Complex.divide(Complex.new(:neg_infinity), -1) == Complex.new(:infinity, 0)

    assert Complex.divide(-1, :infinity) == 0
    assert Complex.divide(0, :infinity) == 0
    assert Complex.divide(1, :infinity) == 0

    assert Complex.divide(-1, :neg_infinity) == 0
    assert Complex.divide(0, :neg_infinity) == 0
    assert Complex.divide(1, :neg_infinity) == 0

    assert Complex.divide(-1, Complex.new(:infinity)) == Complex.new(0)
    assert Complex.divide(0, Complex.new(:infinity)) == Complex.new(0)
    assert Complex.divide(1, Complex.new(:infinity)) == Complex.new(0)

    assert Complex.divide(-1, Complex.new(:neg_infinity)) == Complex.new(0)
    assert Complex.divide(0, Complex.new(:neg_infinity)) == Complex.new(0)
    assert Complex.divide(1, Complex.new(:neg_infinity)) == Complex.new(0)
  end

  test "Math functions" do
    a = Complex.new(1.0, 2.0)
    b = Complex.new(3.0, 4.0)
    c = Complex.new(-1.0, 2.0)
    d = Complex.new(3.0, -4.0)

    assert_close Complex.abs(a), 2.23606797749979
    assert_close Complex.abs(b), 5
    assert_close Complex.abs(c), 2.23606797749979
    assert_close Complex.abs(d), 5
    assert_close Complex.abs_squared(a), 5
    assert_close Complex.abs_squared(b), 25
    assert_close Complex.abs_squared(c), 5
    assert_close Complex.abs_squared(d), 25
    assert Complex.real(a) == 1
    assert Complex.real(b) == 3.0
    assert Complex.real(c) == -1.0
    assert Complex.real(d) == 3.0
    assert Complex.real(:infinity) == :infinity
    assert Complex.real(:neg_infinity) == :neg_infinity
    assert Complex.real(:nan) == :nan
    assert Complex.real(Complex.new(:infinity)) == :infinity
    assert Complex.real(Complex.new(:neg_infinity)) == :neg_infinity
    assert Complex.real(Complex.new(:nan)) == :nan

    assert Complex.imag(a) == 2.0
    assert Complex.imag(b) == 4.0
    assert Complex.imag(c) == 2.0
    assert Complex.imag(d) == -4.0
    assert Complex.imag(:infinity) == 0
    assert Complex.imag(:neg_infinity) == 0
    assert Complex.imag(:nan) == 0
    assert Complex.imag(Complex.new(1, :infinity)) == :infinity
    assert Complex.imag(Complex.new(1, :neg_infinity)) == :neg_infinity
    assert Complex.imag(Complex.new(1, :nan)) == :nan
    assert_close Complex.phase(a), 1.1071487177940904
    assert_close Complex.conjugate(a), %Complex{re: 1.0, im: -2.0}

    assert Complex.conjugate(Complex.new(:infinity, :infinity)) == %Complex{
             re: :infinity,
             im: :neg_infinity
           }

    assert Complex.conjugate(Complex.new(:infinity, :neg_infinity)) == %Complex{
             re: :infinity,
             im: :infinity
           }

    assert Complex.conjugate(Complex.new(:infinity, :nan)) == %Complex{
             re: :infinity,
             im: :nan
           }

    assert Complex.conjugate(Complex.new(:neg_infinity, :infinity)) == %Complex{
             re: :neg_infinity,
             im: :neg_infinity
           }

    assert Complex.conjugate(Complex.new(:neg_infinity, :neg_infinity)) == %Complex{
             re: :neg_infinity,
             im: :infinity
           }

    assert Complex.conjugate(Complex.new(:neg_infinity, :nan)) == %Complex{
             re: :neg_infinity,
             im: :nan
           }

    assert Complex.conjugate(Complex.new(:nan, :infinity)) == %Complex{
             re: :nan,
             im: :neg_infinity
           }

    assert Complex.conjugate(Complex.new(:nan, :neg_infinity)) == %Complex{
             re: :nan,
             im: :infinity
           }

    assert Complex.conjugate(Complex.new(:nan, :nan)) == %Complex{re: :nan, im: :nan}
    assert_close Complex.square(a), Complex.multiply(a, a)
    assert_close Complex.square(b), Complex.multiply(b, b)
    assert_close Complex.square(c), Complex.multiply(c, c)
    assert_close Complex.square(d), Complex.multiply(d, d)
  end

  test "Math functions (no upcasts)" do
    a = 1
    b = 2
    c = 3
    d = 4

    assert_close Complex.abs(a), 1
    assert_close Complex.abs(b), 2
    assert_close Complex.abs(c), 3
    assert_close Complex.abs(d), 4
    assert_close Complex.abs_squared(a), 1
    assert_close Complex.abs_squared(b), 4
    assert_close Complex.abs_squared(c), 9
    assert_close Complex.abs_squared(d), 16
    assert Complex.real(a) == 1
    assert Complex.real(b) == 2
    assert Complex.real(c) == 3
    assert Complex.real(d) == 4
    assert Complex.imag(a) == 0
    assert Complex.imag(b) == 0
    assert Complex.imag(c) == 0
    assert Complex.imag(d) == 0
    assert_close Complex.phase(a), 0
    assert_close Complex.phase(-a), :math.pi()

    assert_close Complex.phase(:infinity), 0
    assert_close Complex.phase(:neg_infinity), :math.pi()
    assert Complex.phase(:nan) == :nan

    for im <- [0, 0.0],
        {re, phase} <- [{:infinity, 0}, {:neg_infinity, :math.pi()}, {:nan, :nan}] do
      assert Complex.phase(Complex.new(re, im)) == phase
    end

    assert Complex.phase(Complex.new(:infinity, :infinity)) == :math.pi() / 4
    assert Complex.phase(Complex.new(:infinity, :nan)) == :nan
    assert Complex.phase(Complex.new(:infinity, :neg_infinity)) == -:math.pi() / 4
    assert Complex.phase(Complex.new(:infinity, -1)) == 0
    assert Complex.phase(Complex.new(:infinity, 1)) == 0
    assert Complex.phase(Complex.new(:neg_infinity, :infinity)) == 3 * :math.pi() / 4
    assert Complex.phase(Complex.new(:neg_infinity, :nan)) == :nan
    assert Complex.phase(Complex.new(:neg_infinity, :neg_infinity)) == 5 * :math.pi() / 4
    assert Complex.phase(Complex.new(:neg_infinity, -1)) == :math.pi()
    assert Complex.phase(Complex.new(:neg_infinity, 1)) == :math.pi()

    assert Complex.phase(Complex.new(1, :infinity)) == :math.pi() / 2
    assert Complex.phase(Complex.new(1, :neg_infinity)) == -:math.pi() / 2
    assert Complex.phase(Complex.new(1, :nan)) == :nan

    assert_close Complex.conjugate(a), a
    assert_close Complex.square(a), a * a
    assert_close Complex.square(d), d * d
  end

  test "abs (non finite)" do
    assert Complex.abs(:nan) == :nan
    assert Complex.abs(:infinity) == :infinity
    assert Complex.abs(:neg_infinity) == :infinity
    assert Complex.abs(Complex.new(:nan, :nan)) == :nan

    assert Complex.abs_squared(:nan) == :nan
    assert Complex.abs_squared(:infinity) == :infinity
    assert Complex.abs_squared(:neg_infinity) == :infinity
    assert Complex.abs_squared(Complex.new(:nan, :nan)) == :nan

    values = [:infinity, :neg_infinity, 0, 1, -1]

    for x <- values, y <- values, is_atom(x) or is_atom(y) do
      assert Complex.abs(Complex.new(x, y)) == :infinity
      assert Complex.abs_squared(Complex.new(x, y)) == :infinity
    end
  end

  test "Exp and logs" do
    a = Complex.new(1.0, 2.0)
    b = Complex.new(3.0, 4.0)
    assert_close Complex.exp(a), %Complex{re: -1.1312043837568135, im: 2.4717266720048188}
    assert_close Complex.log(a), %Complex{re: 0.8047189562170503, im: 1.1071487177940904}
    assert_close Complex.log10(a), %Complex{re: 0.3494850021680094, im: 0.480828578784234}
    assert_close Complex.log2(a), %Complex{re: 1.1609640474436813, im: 1.5972779646881088}
    assert_close Complex.pow(a, b), %Complex{re: 0.129009594074467, im: 0.03392409290517014}
  end

  test "Exp and logs (no upcasts)" do
    a = 3
    b = 4.0
    assert_close Complex.exp(a), :math.exp(a)
    assert_close Complex.log(a), :math.log(a)
    assert_close Complex.log10(a), :math.log10(a)
    assert_close Complex.log2(a), :math.log2(a)
    assert_close Complex.pow(a, b), :math.pow(a, b)

    assert Complex.log(0) == :neg_infinity
    assert Complex.log10(0) == :neg_infinity
    assert Complex.log2(0) == :neg_infinity

    assert Complex.log(:infinity) == :infinity
    assert Complex.log10(:infinity) == :infinity
    assert Complex.log2(:infinity) == :infinity

    assert Complex.log(-1) == :nan
    assert Complex.log10(-1) == :nan
    assert Complex.log2(-1) == :nan

    assert Complex.log(:nan) == :nan
    assert Complex.log10(:nan) == :nan
    assert Complex.log2(:nan) == :nan
  end

  test "pow (non-finite)" do
    assert Complex.pow(:nan, :rand.uniform()) == :nan
    assert Complex.pow(:rand.uniform(), :nan) == :nan

    assert Complex.pow(:infinity, 2) == :infinity
    assert Complex.pow(:infinity, 0) == 1
    assert Complex.pow(:infinity, 0.0) == 1
    assert Complex.pow(:infinity, -2) == 0

    assert Complex.pow(:neg_infinity, 2) == :infinity
    assert Complex.pow(:neg_infinity, 3) == :neg_infinity
    assert Complex.pow(:neg_infinity, 0) == 1
    assert Complex.pow(:neg_infinity, 0.0) == 1
    assert Complex.pow(:neg_infinity, -2) == 0

    assert Complex.pow(10, :infinity) == :infinity
    assert Complex.pow(:infinity, :infinity) == :infinity
    assert Complex.pow(:neg_infinity, :infinity) == :infinity
    assert Complex.pow(10, :neg_infinity) == 0
    assert Complex.pow(:infinity, :neg_infinity) == 0
    assert Complex.pow(:neg_infinity, :neg_infinity) == 0
    assert Complex.pow(0, :neg_infinity) == :infinity
    assert Complex.pow(0, :infinity) == 0
    assert Complex.pow(Complex.new(0, 0), :neg_infinity) == :infinity
    assert Complex.pow(Complex.new(0, 0), :infinity) == 0
  end

  test "log, log10, log2 (non-finite)" do
    assert Complex.log(:infinity) == :infinity
    assert Complex.log(:neg_infinity) == Complex.new(:infinity, :math.pi())
    assert Complex.log(:nan) == :nan

    assert Complex.log(Complex.new(:infinity, :infinity)) ==
             Complex.new(:infinity, :math.pi() / 4)

    assert Complex.log(Complex.new(:infinity, :neg_infinity)) ==
             Complex.new(:infinity, -:math.pi() / 4)

    assert Complex.log(Complex.new(:neg_infinity, :infinity)) ==
             Complex.new(:infinity, 3 * :math.pi() / 4)

    assert Complex.log(Complex.new(:neg_infinity, :neg_infinity)) ==
             Complex.new(:infinity, -3 * :math.pi() / 4)

    assert Complex.log(Complex.new(:infinity, 1)) ==
             Complex.new(:infinity, 0)

    assert Complex.log(Complex.new(:infinity, 0)) ==
             Complex.new(:infinity, 0)

    assert Complex.log(Complex.new(:infinity, -1)) ==
             Complex.new(:infinity, 0)

    assert Complex.log(Complex.new(:neg_infinity, 1)) ==
             Complex.new(:infinity, :math.pi())

    assert Complex.log(Complex.new(:neg_infinity, 0)) ==
             Complex.new(:infinity, :math.pi())

    assert Complex.log(Complex.new(:neg_infinity, -1)) ==
             Complex.new(:infinity, :math.pi())

    assert Complex.log10(:infinity) == :infinity
    assert Complex.log10(:neg_infinity) == Complex.new(:infinity, :math.pi() / :math.log(10))
    assert Complex.log10(:nan) == :nan

    assert Complex.log10(Complex.new(:infinity, :infinity)) ==
             Complex.new(:infinity, :math.pi() / 4 / :math.log(10))

    assert Complex.log10(Complex.new(:infinity, :neg_infinity)) ==
             Complex.new(:infinity, -:math.pi() / 4 / :math.log(10))

    assert Complex.log10(Complex.new(:neg_infinity, :infinity)) ==
             Complex.new(:infinity, 3 * :math.pi() / 4 / :math.log(10))

    assert Complex.log10(Complex.new(:neg_infinity, :neg_infinity)) ==
             Complex.new(:infinity, -3 * :math.pi() / 4 / :math.log(10))

    assert Complex.log10(Complex.new(:infinity, 1)) ==
             Complex.new(:infinity, 0)

    assert Complex.log10(Complex.new(:infinity, 0)) ==
             Complex.new(:infinity, 0)

    assert Complex.log10(Complex.new(:infinity, -1)) ==
             Complex.new(:infinity, 0)

    assert Complex.log10(Complex.new(:neg_infinity, 1)) ==
             Complex.new(:infinity, :math.pi() / :math.log(10))

    assert Complex.log10(Complex.new(:neg_infinity, 0)) ==
             Complex.new(:infinity, :math.pi() / :math.log(10))

    assert Complex.log10(Complex.new(:neg_infinity, -1)) ==
             Complex.new(:infinity, :math.pi() / :math.log(10))

    assert Complex.log2(:infinity) == :infinity
    assert Complex.log2(:neg_infinity) == Complex.new(:infinity, :math.pi() / :math.log(2))
    assert Complex.log2(:nan) == :nan

    assert Complex.log2(Complex.new(:infinity, :infinity)) ==
             Complex.new(:infinity, :math.pi() / 4 / :math.log(2))

    assert Complex.log2(Complex.new(:infinity, :neg_infinity)) ==
             Complex.new(:infinity, -:math.pi() / 4 / :math.log(2))

    assert Complex.log2(Complex.new(:neg_infinity, :infinity)) ==
             Complex.new(:infinity, 3 * :math.pi() / 4 / :math.log(2))

    assert Complex.log2(Complex.new(:neg_infinity, :neg_infinity)) ==
             Complex.new(:infinity, -3 * :math.pi() / 4 / :math.log(2))

    assert Complex.log2(Complex.new(:infinity, 1)) ==
             Complex.new(:infinity, 0)

    assert Complex.log2(Complex.new(:infinity, 0)) ==
             Complex.new(:infinity, 0)

    assert Complex.log2(Complex.new(:infinity, -1)) ==
             Complex.new(:infinity, 0)

    assert Complex.log2(Complex.new(:neg_infinity, 1)) ==
             Complex.new(:infinity, :math.pi() / :math.log(2))

    assert Complex.log2(Complex.new(:neg_infinity, 0)) ==
             Complex.new(:infinity, :math.pi() / :math.log(2))

    assert Complex.log2(Complex.new(:neg_infinity, -1)) ==
             Complex.new(:infinity, :math.pi() / :math.log(2))
  end

  test "atan2 (non-finite)" do
    assert Complex.atan2(:infinity, :infinity) == :math.pi() / 4
    assert Complex.atan2(:infinity, :neg_infinity) == 3 * :math.pi() / 4
    assert Complex.atan2(:infinity, :nan) == :nan
    assert Complex.atan2(:infinity, -1) == :math.pi() / 2
    assert Complex.atan2(:infinity, 0) == :math.pi() / 2
    assert Complex.atan2(:infinity, 1) == :math.pi() / 2
    assert Complex.atan2(:neg_infinity, :infinity) == -:math.pi() / 4
    assert Complex.atan2(:neg_infinity, :neg_infinity) == -3 * :math.pi() / 4
    assert Complex.atan2(:neg_infinity, :nan) == :nan
    assert Complex.atan2(:neg_infinity, -1) == -:math.pi() / 2
    assert Complex.atan2(:neg_infinity, 0) == -:math.pi() / 2
    assert Complex.atan2(:neg_infinity, 1) == -:math.pi() / 2
    assert Complex.atan2(:nan, :infinity) == :nan
    assert Complex.atan2(:nan, :neg_infinity) == :nan
    assert Complex.atan2(:nan, -1) == :nan
    assert Complex.atan2(:nan, 0) == :nan
    assert Complex.atan2(:nan, 1) == :nan
  end

  test "sqrt (non-finite)" do
    assert Complex.sqrt(:infinity) == :infinity
    assert Complex.sqrt(:neg_infinity) == :nan
    assert Complex.sqrt(:nan) == :nan

    assert Complex.sqrt(Complex.new(:infinity, 1)) == Complex.new(:infinity)
    assert Complex.sqrt(Complex.new(:neg_infinity, 1)) == Complex.new(0, :infinity)
    assert Complex.sqrt(Complex.new(:infinity, :infinity)) == Complex.new(:infinity, :infinity)

    assert Complex.sqrt(Complex.new(:neg_infinity, :infinity)) ==
             Complex.new(:infinity, :infinity)

    assert Complex.sqrt(Complex.new(:infinity, :neg_infinity)) ==
             Complex.new(:infinity, :neg_infinity)

    assert Complex.sqrt(Complex.new(:neg_infinity, :neg_infinity)) ==
             Complex.new(:infinity, :neg_infinity)
  end

  test "exp (non-finite)" do
    assert Complex.exp(:infinity) == :infinity
    assert Complex.exp(:neg_infinity) == 0
    assert Complex.exp(:nan) == :nan

    for x <- [:infinity, :neg_infinity, :nan, -1, 0, 1] do
      assert Complex.exp(Complex.new(:neg_infinity, x)) == 0
    end

    assert Complex.exp(Complex.new(:infinity, :nan)) == Complex.new(:infinity, :nan)
    assert Complex.exp(Complex.new(:infinity, -1)) == Complex.new(:infinity, :neg_infinity)
    assert Complex.exp(Complex.new(:infinity, 0)) == Complex.new(:infinity, :nan)
    assert Complex.exp(Complex.new(:infinity, 1)) == Complex.new(:infinity, :infinity)

    for x <- [:infinity, :neg_infinity, :nan] do
      assert Complex.exp(Complex.new(0, x)) == Complex.new(:nan, :nan)
    end

    assert Complex.exp(Complex.new(:nan, 1)) == Complex.new(:nan, :nan)
  end

  test "Trig functions" do
    a = Complex.new(1.0, 2.0)
    assert_close Complex.sin(a), %Complex{re: 3.165778513216168, im: 1.9596010414216063}
    assert_close Complex.asin(a), %Complex{re: 0.4270785863924759, im: 1.5285709194809978}
    assert_close Complex.cos(a), %Complex{re: 2.0327230070196656, im: -3.0518977991518}
    assert_close Complex.acos(a), %Complex{re: 1.1437177404024206, im: -1.528570919480998}
    assert_close Complex.tan(a), %Complex{re: 0.03381282607989661, im: 1.0147936161466335}
    assert_close Complex.atan(a), %Complex{re: 1.3389725222944935, im: 0.4023594781085251}
    assert_close Complex.cot(a), %Complex{re: 0.032797755533752464, im: -0.9843292264581909}
    assert_close Complex.acot(a), %Complex{re: 0.23182380450040307, im: -0.4023594781085251}
    assert_close Complex.sec(a), %Complex{re: 0.1511762982655772, im: 0.22697367539372157}
    assert_close Complex.asec(a), %Complex{re: 1.384478272687081, im: 0.39656823011232895}
    assert_close Complex.csc(a), %Complex{re: 0.22837506559968654, im: -0.1413630216124078}
    assert_close Complex.acsc(a), %Complex{re: 0.18631805410781554, im: -0.3965682301123289}
  end

  test "Trig functions (no upcasts)" do
    a = 3
    assert_close Complex.sin(a), 0.141120
    assert_close Complex.asin(0.141120), :math.pi() - a
    assert_close Complex.cos(a), -0.989992
    assert_close Complex.acos(-0.989992), a
    assert_close Complex.tan(a), -0.1425466
    assert_close Complex.atan(-0.1425466), a - :math.pi()
    assert_close Complex.cot(a), -7.015252
    assert_close Complex.acot(-7.015252), a - :math.pi()
    assert_close Complex.sec(a), -1.010108
    assert_close Complex.asec(-1.010108), a
    assert_close Complex.csc(a), 7.086167
    assert_close Complex.acsc(7.086167), :math.pi() - a
  end

  test "sin (non-finite)" do
    assert Complex.sin(Complex.new(:nan)) == Complex.new(:nan, :nan)
    assert Complex.sin(Complex.new(:infinity)) == Complex.new(:nan, :nan)
    assert Complex.sin(Complex.new(:neg_infinity)) == Complex.new(:nan, :nan)
    assert Complex.sin(Complex.new(0, :nan)) == Complex.new(:nan, :nan)
    assert Complex.sin(Complex.new(0, :infinity)) == Complex.new(:nan, :infinity)
    assert Complex.sin(Complex.new(0, :neg_infinity)) == Complex.new(:nan, :neg_infinity)

    assert Complex.sin(Complex.new(-1, :infinity)) == Complex.new(:neg_infinity, :infinity)
    assert Complex.sin(Complex.new(1, :infinity)) == Complex.new(:infinity, :infinity)
    assert Complex.sin(Complex.new(-1, :neg_infinity)) == Complex.new(:infinity, :neg_infinity)
    assert Complex.sin(Complex.new(1, :neg_infinity)) == Complex.new(:neg_infinity, :neg_infinity)
  end

  test "cos (non-finite)" do
    assert Complex.cos(Complex.new(:nan)) == Complex.new(:nan, :nan)
    assert Complex.cos(Complex.new(:infinity)) == Complex.new(:nan, :nan)
    assert Complex.cos(Complex.new(:neg_infinity)) == Complex.new(:nan, :nan)
    assert Complex.cos(Complex.new(0, :nan)) == Complex.new(:nan, :nan)
    assert Complex.cos(Complex.new(0, :infinity)) == Complex.new(:infinity, :nan)
    assert Complex.cos(Complex.new(0, :neg_infinity)) == Complex.new(:infinity, :nan)

    assert Complex.cos(Complex.new(-1, :infinity)) == Complex.new(:infinity, :infinity)
    assert Complex.cos(Complex.new(1, :infinity)) == Complex.new(:infinity, :neg_infinity)
    assert Complex.cos(Complex.new(-1, :neg_infinity)) == Complex.new(:infinity, :neg_infinity)
    assert Complex.cos(Complex.new(1, :neg_infinity)) == Complex.new(:infinity, :infinity)
  end

  test "Hyperbolic functions" do
    a = Complex.new(1.0, 2.0)
    assert_close Complex.sinh(a), %Complex{re: -0.48905625904129363, im: 1.4031192506220405}
    assert_close Complex.asinh(a), %Complex{re: 1.4693517443681852, im: 1.0634400235777521}
    assert_close Complex.cosh(a), %Complex{re: -0.64214812471552, im: 1.0686074213827783}
    assert_close Complex.acosh(a), %Complex{re: 1.528570919480998, im: 1.1437177404024206}
    assert_close Complex.tanh(a), %Complex{re: 1.16673625724092, im: -0.24345820118572528}
    assert_close Complex.atanh(a), %Complex{re: 0.17328679513998635, im: 1.1780972450961724}
    assert_close Complex.sech(a), %Complex{re: -0.41314934426694006, im: -0.687527438655479}
    assert_close Complex.asech(a), %Complex{re: 0.3965682301123289, im: -1.3844782726870812}
    assert_close Complex.csch(a), %Complex{re: -0.22150093085050943, im: -0.6354937992538999}
    assert_close Complex.acsch(a), %Complex{re: 0.2156124185558298, im: -0.4015863916678061}
    assert_close Complex.coth(a), %Complex{re: 0.8213297974938518, im: 0.171383612909185}
    assert_close Complex.acoth(a), %Complex{re: 0.1732867951399863, im: -0.39269908169872414}
  end

  test "Hyperbolic functions (no upcasts)" do
    a = 3
    assert_close Complex.sinh(a), 10.017874
    assert_close Complex.asinh(10.017874), a
    assert_close Complex.cosh(a), 10.06766
    assert_close Complex.acosh(10.06766), a
    assert_close Complex.tanh(a), 0.9950547
    assert_close Complex.atanh(0.9950547), a
    assert_close Complex.sech(a), 0.099327
    assert_close Complex.asech(0.099327), a
    assert_close Complex.csch(a), 0.0998215
    assert_close Complex.acsch(0.0998215), a
    assert_close Complex.coth(a), 1.0049698
    assert_close Complex.acoth(1.0049698), a
  end

  test "trig and hyperbolic trig functions non-finite" do
    for n <- [:infinity, :neg_infinity, :nan] do
      assert Complex.asin(n) == :nan
      assert Complex.asin(Complex.new(n)) == Complex.new(:nan, :nan)
      assert Complex.asin(Complex.new(n, 1)) == Complex.new(:nan, :nan)
      assert Complex.asin(Complex.new(1, n)) == Complex.new(:nan, :nan)
      assert Complex.asin(Complex.new(0, n)) == Complex.new(:nan, :nan)
      assert Complex.acos(n) == :nan
      assert Complex.acos(Complex.new(n)) == Complex.new(:nan, :nan)
      assert Complex.acos(Complex.new(n, 1)) == Complex.new(:nan, :nan)
      assert Complex.acos(Complex.new(1, n)) == Complex.new(:nan, :nan)
      assert Complex.acos(Complex.new(0, n)) == Complex.new(:nan, :nan)
      assert Complex.tan(Complex.new(n)) == Complex.new(:nan, :nan)
      assert Complex.tan(Complex.new(n, 1)) == Complex.new(:nan, :nan)
      assert Complex.tan(Complex.new(1, n)) == Complex.new(:nan, :nan)
      assert Complex.tan(Complex.new(0, n)) == Complex.new(:nan, :nan)
      assert Complex.cot(Complex.new(n)) == Complex.new(:nan, :nan)
      assert Complex.cot(Complex.new(n, 1)) == Complex.new(:nan, :nan)
      assert Complex.cot(Complex.new(1, n)) == Complex.new(:nan, :nan)
      assert Complex.cot(Complex.new(0, n)) == Complex.new(:nan, :nan)
      assert Complex.asinh(n) == n
    end

    assert Complex.acot(:infinity) == 0
    assert Complex.acot(:neg_infinity) == :math.pi()
    assert Complex.acot(:nan) == :nan

    assert Complex.asec(:infinity) == :math.pi() / 2
    assert Complex.asec(:neg_infinity) == 3 * :math.pi() / 2
    assert Complex.asec(:nan) == :nan

    assert Complex.acsc(:infinity) == 0
    assert Complex.acsc(:neg_infinity) == -:math.pi()
    assert Complex.acsc(:nan) == :nan

    assert Complex.acosh(:infinity) == :infinity
    assert Complex.acosh(:nan) == :nan

    assert_raise ArithmeticError, "Complex.acosh(:neg_infinity) is undefined", fn ->
      Complex.acosh(:neg_infinity)
    end
  end

  for {m, f} <- [{Kernel, :inspect}, {String.Chars.Complex, :to_string}, {Complex, :to_string}] do
    test "#{m}.#{f}/1" do
      assert "1.0+1.0i" == apply(unquote(m), unquote(f), [Complex.new(1.0, 1.0)])
      assert "1.0-1.0i" == apply(unquote(m), unquote(f), [Complex.new(1.0, -1.0)])
      assert "1.0+0.0i" == apply(unquote(m), unquote(f), [Complex.new(1.0, 0.0)])
      assert "1.0+0.0i" == apply(unquote(m), unquote(f), [Complex.new(1.0, -0.0)])

      assert "Inf+0.0i" == apply(unquote(m), unquote(f), [Complex.new(:infinity, -0.0)])
      assert "-Inf+0.0i" == apply(unquote(m), unquote(f), [Complex.new(:neg_infinity, -0.0)])
      assert "NaN+0.0i" == apply(unquote(m), unquote(f), [Complex.new(:nan, -0.0)])
      assert "0.0+Infi" == apply(unquote(m), unquote(f), [Complex.new(0.0, :infinity)])
      assert "0.0-Infi" == apply(unquote(m), unquote(f), [Complex.new(0.0, :neg_infinity)])
      assert "0.0+NaNi" == apply(unquote(m), unquote(f), [Complex.new(0.0, :nan)])

      assert "1.0-1.0i" ==
               apply(unquote(m), unquote(f), [Complex.new(1.0, 1.0) |> Complex.conjugate()])

      assert "1.0+1.0i" ==
               apply(unquote(m), unquote(f), [Complex.new(1.0, -1.0) |> Complex.conjugate()])

      assert "1.0+0.0i" ==
               apply(unquote(m), unquote(f), [Complex.new(1.0, 0.0) |> Complex.conjugate()])

      assert "1.0+0.0i" ==
               apply(unquote(m), unquote(f), [Complex.new(1.0, -0.0) |> Complex.conjugate()])
    end
  end

  test "erf" do
    assert_close Complex.erf(-1), -0.8427
    assert Complex.erf(0) == 0
    assert_close Complex.erf(1), 0.8427
    assert Complex.erf(:nan) == :nan
    assert Complex.erf(:infinity) == 1
    assert Complex.erf(:neg_infinity) == -1
  end

  test "erfc" do
    assert_close Complex.erfc(-1), 1.8427
    assert Complex.erfc(0) == 1
    assert_close Complex.erfc(1), 0.15730
    assert Complex.erfc(:nan) == :nan
    assert Complex.erfc(:infinity) == 0
    assert Complex.erfc(:neg_infinity) == 2
  end

  test "erf_inv" do
    assert_close(Complex.erf_inv(-0.8427007929497149), -1)
    assert Complex.erf_inv(0) == 0
    assert_close(Complex.erf_inv(0.8427007929497149), 1)
    assert Complex.erf_inv(:nan) == :nan
    assert Complex.erf_inv(1) == :infinity
    assert Complex.erf_inv(-1) == :neg_infinity
  end

  test "complex casting - add/2" do
    assert Complex.add(:infinity, 1) == :infinity
    assert Complex.add(:neg_infinity, 1) == :neg_infinity

    assert Complex.add(:infinity, Complex.new(:neg_infinity, 0)) == Complex.new(:nan, 0)
    assert Complex.add(:neg_infinity, Complex.new(:infinity, 0)) == Complex.new(:nan, 0)
    assert Complex.add(:infinity, Complex.new(:infinity, 0)) == Complex.new(:infinity, 0)

    assert Complex.add(:neg_infinity, Complex.new(:neg_infinity, 0)) ==
             Complex.new(:neg_infinity, 0)

    assert Complex.add(Complex.new(:neg_infinity, 0), :infinity) == Complex.new(:nan, 0)
    assert Complex.add(Complex.new(:infinity, 0), :neg_infinity) == Complex.new(:nan, 0)
    assert Complex.add(Complex.new(:infinity, 0), :infinity) == Complex.new(:infinity, 0)

    assert Complex.add(Complex.new(:neg_infinity, 0), :neg_infinity) ==
             Complex.new(:neg_infinity, 0)

    assert Complex.add(:infinity, Complex.new(1, 2)) == Complex.new(:infinity, 2)
    assert Complex.add(Complex.new(1, 2), :infinity) == Complex.new(:infinity, 2)
  end

  test "complex casting - subtract/2" do
    assert Complex.subtract(:infinity, 1) == :infinity
    assert Complex.subtract(:neg_infinity, 1) == :neg_infinity

    assert Complex.subtract(:infinity, Complex.new(:neg_infinity, 0)) == Complex.new(:infinity, 0)

    assert Complex.subtract(:neg_infinity, Complex.new(:infinity, 0)) ==
             Complex.new(:neg_infinity, 0)

    assert Complex.subtract(:infinity, Complex.new(:infinity, 0)) == Complex.new(:nan, 0)

    assert Complex.subtract(:neg_infinity, Complex.new(:neg_infinity, 0)) ==
             Complex.new(:nan, 0)

    assert Complex.subtract(Complex.new(:neg_infinity, 0), :infinity) ==
             Complex.new(:neg_infinity, 0)

    assert Complex.subtract(Complex.new(:infinity, 0), :neg_infinity) == Complex.new(:infinity, 0)
    assert Complex.subtract(Complex.new(:infinity, 0), :infinity) == Complex.new(:nan, 0)

    assert Complex.subtract(Complex.new(:neg_infinity, 0), :neg_infinity) ==
             Complex.new(:nan, 0)

    assert Complex.subtract(:infinity, Complex.new(1, 2)) == Complex.new(:infinity, -2)
    assert Complex.subtract(Complex.new(1, 2), :infinity) == Complex.new(:neg_infinity, 2)
  end

  test "complex casting - multiply/2" do
    assert Complex.multiply(1, 1) == 1
    assert Complex.multiply(:infinity, Complex.new(1, 1)) == Complex.new(:infinity, :infinity)

    assert Complex.multiply(:neg_infinity, Complex.new(1, 1)) ==
             Complex.new(:neg_infinity, :neg_infinity)

    assert Complex.multiply(:infinity, :infinity) == :infinity
    assert Complex.multiply(:neg_infinity, :infinity) == :neg_infinity
    assert Complex.multiply(:infinity, :neg_infinity) == :neg_infinity
    assert Complex.multiply(:neg_infinity, :neg_infinity) == :infinity

    # First direction
    assert Complex.multiply(Complex.new(:infinity, 1), :infinity) == Complex.new(:infinity, :nan)

    assert Complex.multiply(Complex.new(:neg_infinity, 1), :neg_infinity) ==
             Complex.new(:infinity, :nan)

    assert Complex.multiply(Complex.new(:infinity, 1), :neg_infinity) ==
             Complex.new(:neg_infinity, :nan)

    assert Complex.multiply(Complex.new(:neg_infinity, 1), :infinity) ==
             Complex.new(:neg_infinity, :nan)

    assert Complex.multiply(Complex.new(1, :infinity), :infinity) == Complex.new(:nan, :infinity)

    assert Complex.multiply(Complex.new(1, :infinity), :neg_infinity) ==
             Complex.new(:nan, :neg_infinity)

    assert Complex.multiply(Complex.new(1, :neg_infinity), :neg_infinity) ==
             Complex.new(:nan, :infinity)

    assert Complex.multiply(Complex.new(1, :neg_infinity), :infinity) ==
             Complex.new(:nan, :neg_infinity)

    # Second direction
    assert Complex.multiply(:infinity, Complex.new(:infinity, 1)) == Complex.new(:infinity, :nan)

    assert Complex.multiply(:neg_infinity, Complex.new(:neg_infinity, 1)) ==
             Complex.new(:infinity, :nan)

    assert Complex.multiply(:neg_infinity, Complex.new(:infinity, 1)) ==
             Complex.new(:neg_infinity, :nan)

    assert Complex.multiply(:infinity, Complex.new(:neg_infinity, 1)) ==
             Complex.new(:neg_infinity, :nan)

    assert Complex.multiply(:infinity, Complex.new(1, :infinity)) == Complex.new(:nan, :infinity)

    assert Complex.multiply(:neg_infinity, Complex.new(1, :infinity)) ==
             Complex.new(:nan, :neg_infinity)

    assert Complex.multiply(:neg_infinity, Complex.new(1, :neg_infinity)) ==
             Complex.new(:nan, :infinity)

    assert Complex.multiply(:infinity, Complex.new(1, :neg_infinity)) ==
             Complex.new(:nan, :neg_infinity)

    # NaNs
    assert Complex.multiply(:nan, 1) == :nan
    assert Complex.multiply(1, :nan) == :nan

    ## First direction
    assert Complex.multiply(:nan, Complex.new(1, 1)) == Complex.new(:nan, :nan)
    assert Complex.multiply(1, Complex.new(:nan, 1)) == Complex.new(:nan, :nan)
    assert Complex.multiply(1, Complex.new(1, :nan)) == Complex.new(:nan, :nan)

    ## Second direction
    assert Complex.multiply(Complex.new(1, 1), :nan) == Complex.new(:nan, :nan)
    assert Complex.multiply(Complex.new(:nan, 1), 1) == Complex.new(:nan, :nan)
    assert Complex.multiply(Complex.new(1, :nan), 1) == Complex.new(:nan, :nan)
  end

  test "complex casting - divide/2" do
    assert Complex.divide(1, 1) == 1

    assert Complex.divide(:infinity, Complex.new(:infinity, 1)) == Complex.new(:nan, :nan)
    assert Complex.divide(:neg_infinity, Complex.new(:neg_infinity, 1)) == Complex.new(:nan, :nan)

    assert Complex.divide(:infinity, Complex.new(1, 1)) == Complex.new(:infinity, :neg_infinity)

    assert Complex.divide(:neg_infinity, Complex.new(1, 1)) ==
             Complex.new(:neg_infinity, :infinity)

    assert Complex.divide(Complex.new(1, 1), :infinity) == Complex.new(0, 0)
    assert Complex.divide(Complex.new(1, 1), :neg_infinity) == Complex.new(0, 0)

    assert Complex.divide(Complex.new(:infinity, 1), :infinity) == Complex.new(:nan, :nan)
    assert Complex.divide(Complex.new(1, :infinity), :infinity) == Complex.new(:nan, :nan)
    assert Complex.divide(Complex.new(:infinity, 1), :neg_infinity) == Complex.new(:nan, :nan)
    assert Complex.divide(Complex.new(1, :infinity), :neg_infinity) == Complex.new(:nan, :nan)

    # NaNs
    assert Complex.divide(:nan, 1) == :nan
    assert Complex.divide(1, :nan) == :nan

    ## First direction
    assert Complex.divide(:nan, Complex.new(1, 1)) == Complex.new(:nan, :nan)
    assert Complex.divide(1, Complex.new(:nan, 1)) == Complex.new(:nan, :nan)
    assert Complex.divide(1, Complex.new(1, :nan)) == Complex.new(:nan, :nan)

    ## Second direction
    assert Complex.divide(Complex.new(1, 1), :nan) == Complex.new(:nan, :nan)
    assert Complex.divide(Complex.new(:nan, 1), 1) == Complex.new(:nan, :nan)
    assert Complex.divide(Complex.new(1, :nan), 1) == Complex.new(:nan, :nan)
  end

  test "complex casting - pow/2"

  test "complex casting - atan2/2"

  defp assert_close(left, right, opts \\ []) do
    eps = opts[:eps] || 1.0e-5

    left
    |> Complex.subtract(right)
    |> Complex.abs()
    |> case do
      val when val <= eps ->
        true

      _ ->
        flunk("""
        Expected the two values to be close.
        left: #{inspect(left)}
        right: #{inspect(right)}
        """)
    end
  end
end
