defmodule Complex.SigilsTest do
  use ExUnit.Case, async: true

  import Complex, only: :sigils

  test "~COMPLEX builds a complex number from integers" do
    assert ~COMPLEX(3 4) == Complex.new(3, 4)
  end

  test "~COMPLEX keeps integers and floats as written" do
    assert ~COMPLEX(3 -4.5) === Complex.new(3, -4.5)
  end

  test "extra whitespace is ignored" do
    assert ~COMPLEX(  3    4  ) == Complex.new(3, 4)
  end

  test "~POLAR builds a complex number from magnitude and angle" do
    assert ~POLAR(2 1.57) == Complex.from_polar(2, 1.57)
  end

  test "it raises on an invalid number" do
    assert_raise ArgumentError, ~r/invalid number/, fn ->
      Complex.sigil_COMPLEX("3 x", [])
    end

    assert_raise ArgumentError, ~r/invalid number/, fn ->
      Complex.sigil_POLAR("3 x", [])
    end
  end

  test "it raises when there are not exactly two numbers" do
    assert_raise MatchError, fn -> Complex.sigil_COMPLEX("3", []) end
    assert_raise MatchError, fn -> Complex.sigil_POLAR("1 2 3", []) end
  end

  test "it raises on modifiers" do
    # apply/3 keeps the compiler's type checker from flagging the call,
    # which can never match the `[]` modifiers clause
    assert_raise FunctionClauseError, fn -> apply(Complex, :sigil_COMPLEX, ["3 4", ~c"x"]) end
    assert_raise FunctionClauseError, fn -> apply(Complex, :sigil_POLAR, ["3 4", ~c"p"]) end
  end
end
