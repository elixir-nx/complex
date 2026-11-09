defmodule Complex.Sigils do
  @moduledoc """
  Sigils to build complex numbers.
  """

  @doc ~S"""
  Handles the sigil `~z` for complex numbers.

  The two numbers are separated by whitespace, mirroring Common Lisp's
  `#c(re im)` reader syntax. Without modifiers they are the real and
  imaginary parts. With the `p` modifier (an extension specific to this
  library) they are the magnitude and the angle (radians).

  ## Examples

      import Complex.Sigils

      ~z(3 4)
      #=> %Complex{re: 3, im: 4}

      ~z(2 1.57)p
      #=> Complex.from_polar(2, 1.57)
  """
  def sigil_z(string, []) do
    [re, im] = parse(string)
    Complex.new(re, im)
  end

  def sigil_z(string, [?p]) do
    [r, phi] = parse(string)
    Complex.from_polar(r, phi)
  end

  defp parse(string), do: string |> String.split() |> Enum.map(&to_number/1)

  defp to_number(str) do
    with :error <- parse_exact(Integer.parse(str)),
         :error <- parse_exact(Float.parse(str)) do
      raise ArgumentError, "invalid number in ~z sigil: #{inspect(str)}"
    end
  end

  defp parse_exact({n, ""}), do: n
  defp parse_exact(_), do: :error
end
