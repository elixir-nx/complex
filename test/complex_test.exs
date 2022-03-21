defmodule ComplexTest do
  use ExUnit.Case
  doctest Complex

  test "Creation" do
    assert Complex.new(1.0, 1.0) == %Complex{re: 1.0, im: 1.0}
    assert Complex.new(1.0, 2.0) == %Complex{re: 1.0, im: 2.0}
    assert Complex.new(2.0) == %Complex{re: 2.0, im: 0.0}

    a = Complex.new(1.0, 2.0)
    assert Complex.getPolar(a) == {2.23606797749979, 1.1071487177940904}
    assert Complex.fromPolar(1.0, :math.pi() / 2) == %Complex{re: 6.123233995736766e-17, im: 1.0}
  end

  test "Parse from string" do
    assert Complex.parse("1.0+1.0i") == %Complex{re: 1.0, im: 1.0}
    # Test various spacing
    assert Complex.parse("1.0 +1.0i") == %Complex{re: 1.0, im: 1.0}
    assert Complex.parse("1.0+ 1.0i") == %Complex{re: 1.0, im: 1.0}
    assert Complex.parse("1.0 + 1.0i") == %Complex{re: 1.0, im: 1.0}

    # We don't handle integer values for re and im
    assert_raise MatchError, fn -> Complex.parse("1+1i") end
    assert_raise MatchError, fn -> Complex.parse("1+1.0i") end
    assert_raise MatchError, fn -> Complex.parse("1.0+1i") end
  end

  test "Arithmetic" do
    a = Complex.new(1.0, 2.0)
    b = Complex.new(3.0, 4.0)
    assert Complex.add(a, b) == %Complex{re: 4.0, im: 6.0}
    assert Complex.sub(a, b) == %Complex{re: -2.0, im: -2.0}
    assert Complex.mult(a, b) == %Complex{re: -5.0, im: 10.0}
    assert Complex.div(a, b) == %Complex{re: 0.44, im: 0.08}
    assert Complex.div(b, a) == %Complex{re: 2.2, im: -0.4}
  end

  test "Math functions" do
    a = Complex.new(1.0, 2.0)
    b = Complex.new(3.0, 4.0)
    c = Complex.new(-1.0, 2.0)
    d = Complex.new(3.0, -4.0)

    assert Complex.abs(a) == 2.23606797749979
    assert Complex.abs(b) == 5
    assert Complex.abs(c) == 2.23606797749979
    assert Complex.abs(d) == 5
    assert Complex.abs_squared(a) == 5
    assert Complex.abs_squared(b) == 25
    assert Complex.abs_squared(c) == 5
    assert Complex.abs_squared(d) == 25
    assert Complex.phase(a) == 1.1071487177940904
    assert Complex.conjugate(a) == %Complex{re: 1.0, im: -2.0}
    assert Complex.square(a) == Complex.mult(a, a)
    assert Complex.square(b) == Complex.mult(b, b)
    assert Complex.square(c) == Complex.mult(c, c)
    assert Complex.square(d) == Complex.mult(d, d)
  end

  test "Exp and logs" do
    a = Complex.new(1.0, 2.0)
    b = Complex.new(3.0, 4.0)
    assert Complex.exp(a) == %Complex{re: -1.1312043837568135, im: 2.4717266720048188}
    assert Complex.ln(a) == %Complex{re: 0.8047189562170503, im: 1.1071487177940904}
    assert Complex.log10(a) == %Complex{re: 0.3494850021680094, im: 0.480828578784234}
    assert Complex.log2(a) == %Complex{re: 1.1609640474436813, im: 1.5972779646881088}
    assert Complex.pow(a, b) == %Complex{re: 0.129009594074467, im: 0.03392409290517014}
  end

  test "Trig functions" do
    a = Complex.new(1.0, 2.0)
    assert Complex.sin(a) == %Complex{re: 3.165778513216168, im: 1.9596010414216063}
    assert Complex.asin(a) == %Complex{re: 0.4270785863924759, im: 1.5285709194809978}
    assert Complex.cos(a) == %Complex{re: 2.0327230070196656, im: -3.0518977991518}
    assert Complex.acos(a) == %Complex{re: 1.1437177404024206, im: -1.528570919480998}
    assert Complex.tan(a) == %Complex{re: 0.03381282607989661, im: 1.0147936161466335}
    assert Complex.atan(a) == %Complex{re: 1.3389725222944935, im: 0.4023594781085251}
    assert Complex.cot(a) == %Complex{re: 0.032797755533752464, im: -0.9843292264581909}
    assert Complex.acot(a) == %Complex{re: 0.23182380450040307, im: -0.4023594781085251}
    assert Complex.sec(a) == %Complex{re: 0.1511762982655772, im: 0.22697367539372157}
    assert Complex.asec(a) == %Complex{re: 1.384478272687081, im: 0.39656823011232895}
    assert Complex.csc(a) == %Complex{re: 0.22837506559968654, im: -0.1413630216124078}
    assert Complex.acsc(a) == %Complex{re: 0.18631805410781554, im: -0.3965682301123289}
  end

  test "Hyperbolic functions" do
    a = Complex.new(1.0, 2.0)
    assert Complex.sinh(a) == %Complex{re: -0.48905625904129363, im: 1.4031192506220405}
    assert Complex.asinh(a) == %Complex{re: 1.4693517443681852, im: 1.0634400235777521}
    assert Complex.cosh(a) == %Complex{re: -0.64214812471552, im: 1.0686074213827783}
    assert Complex.acosh(a) == %Complex{re: 1.528570919480998, im: 1.1437177404024206}
    assert Complex.tanh(a) == %Complex{re: 1.16673625724092, im: -0.24345820118572528}
    assert Complex.atanh(a) == %Complex{re: 0.17328679513998635, im: 1.1780972450961724}
    assert Complex.sech(a) == %Complex{re: -0.41314934426694006, im: -0.687527438655479}
    assert Complex.asech(a) == %Complex{re: 0.3965682301123289, im: -1.3844782726870812}
    assert Complex.csch(a) == %Complex{re: -0.22150093085050943, im: -0.6354937992538999}
    assert Complex.acsch(a) == %Complex{re: 0.2156124185558298, im: -0.4015863916678061}
    assert Complex.coth(a) == %Complex{re: 0.8213297974938518, im: 0.171383612909185}
    assert Complex.acoth(a) == %Complex{re: 0.1732867951399863, im: -0.39269908169872414}
  end
end
