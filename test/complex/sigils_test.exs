defmodule Complex.SigilsTest do
  use ExUnit.Case, async: true

  import Complex.Sigils

  test "cartesian form with integers" do
    assert ~z(3 4) == Complex.new(3, 4)
  end

  test "cartesian form keeps integers and floats as written" do
    assert ~z(3 -4.5) === Complex.new(3, -4.5)
  end

  test "extra whitespace is ignored" do
    assert ~z(  3    4  ) == Complex.new(3, 4)
  end

  test "polar form with the p modifier" do
    assert ~z(2 1.57)p == Complex.from_polar(2, 1.57)
  end

  test "it raises on an invalid number" do
    assert_raise ArgumentError, ~r/invalid number/, fn -> sigil_z("3 x", []) end
  end

  test "it raises when there are not exactly two numbers" do
    assert_raise MatchError, fn -> sigil_z("3", []) end
    assert_raise MatchError, fn -> sigil_z("1 2 3", []) end
  end

  test "it raises on an unknown modifier" do
    assert_raise FunctionClauseError, fn -> sigil_z("3 4", ~c"x") end
  end
end
