defmodule Complex.KernelSigilTest do
  use ExUnit.Case, async: true
  use Complex.Kernel

  test "use Complex.Kernel imports the ~z sigil" do
    assert ~z(3 4) == Complex.new(3, 4)
    assert ~z(2 1.57)p == Complex.from_polar(2, 1.57)
  end
end
