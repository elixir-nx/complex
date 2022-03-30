defmodule Complex.Kernel do
  @moduledoc """
  Provides operator overloading for Elixir operators.
  """

  defmacro __using__(_) do
    has_pow_operator = function_exported?(Kernel, :**, 2)

    operators = [+: 2, -: 2, /: 2, *: 2]

    operators =
      if has_pow_operator do
        [{:**, 2} | operators]
      else
        operators
      end

    quote do
      import Kernel, except: unquote(operators)
      import Complex.Kernel, only: unquote(operators)
    end
  end

  def a + b, do: Complex.add(a, b)
  def a - b, do: Complex.subtract(a, b)
  def a * b, do: Complex.multiply(a, b)
  def a / b, do: Complex.divide(a, b)

  def unquote(:**)(a, b), do: Complex.power(a, b)
end
