defmodule Complex.CompileTimeChecks do
  @moduledoc false

  def math_fun_supported?(fun, arity) do
    Code.ensure_loaded?(:math) and check_math_fun_supported?(fun, arity)
  end

  defp check_math_fun_supported?(fun, arity) do
    args =
      case {fun, arity} do
        {:atan, 1} -> [3.14]
        {:atanh, 1} -> [0.9]
        {_, 1} -> [1.0]
        {_, 2} -> [1.0, 1.0]
      end

    _ = apply(:math, fun, args)
    true
  rescue
    UndefinedFunctionError ->
      false
  end
end
