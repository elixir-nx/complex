defmodule Complex do
  @moduledoc """
  *Complex* is a library for types and mathematical functions for complex
  numbers.

  Each complex number is represented as a structure holding the real and
  imaginary part.  There are functions for creation and manipulation of
  them.  Unfortunately since there is no operator overloading in Elixir the
  math functions (add, subtract, etc.) are implemented as add/2, subtract/2, etc.

  ## Examples
      iex> Complex.new(3, 4)
      %Complex{im: 4, re: 3}

      iex> Complex.imag()
      %Complex{im: 1.0, re: 0.0}
  """

  @vsn 2

  import Kernel, except: [abs: 1]
  import Complex.CompileTimeChecks

  @typedoc """
  General type for complex numbers
  """
  @type complex :: %Complex{re: number, im: number}
  defstruct re: 0, im: 0

  defimpl Inspect do
    def inspect(val, opts),
      do: Complex.to_string(val, opts.custom_options[:imaginary_constant] || "i")
  end

  def to_string(%Complex{re: re, im: im}, imaginary_constant \\ "i") do
    if im < 0 do
      "#{re} - #{abs(im)}#{imaginary_constant}"
    else
      "#{re} + #{im}#{imaginary_constant}"
    end
  end

  @doc """
  Returns a new complex with specified real and imaginary components.  The
  imaginary part defaults to zero so a "real" number can be created with new/1

  #### See also
  [imag/0](#imag/0) [from_polar/2](#from_polar/2)

  #### Examples
      iex> Complex.new(3, 4)
      %Complex{im: 4, re: 3}

      iex> Complex.new(2)
      %Complex{im: 0, re: 2}
  """
  @spec new(number, number) :: complex
  def new(re, im \\ 0), do: %Complex{re: re, im: im}

  @doc """
  Parses a complex number from a string.  The values of the real and imaginary
  parts must be represented by a float, including decimal and at least one
  trailing digit (e.g. 1.2, 0.4).

  #### See also
  [new/2](#new/2)

  #### Examples
      iex> Complex.parse("1.1+2.2i")
      %Complex{im: 2.2, re: 1.1}
  """
  @spec parse(String.t()) :: complex
  def parse(str) do
    [_, real, imag] = Regex.run(~r/([-]?\d+\.\d+)\s*\+\s*([-]?\d+\.\d+)i/, str)

    Complex.new(String.to_float(real), String.to_float(imag))
  end

  @doc """
  Returns a new complex representing the pure imaginary number sqrt(-1).

  #### See also
  [new/2](#new/2) [from_polar/2](#from_polar/2)

  #### Examples
      iex> Complex.imag()
      %Complex{im: 1.0, re: 0.0}
  """
  @spec imag() :: complex
  def imag(), do: %Complex{re: 0.0, im: 1.0}

  @doc """
  Returns a new complex number described by the supplied polar coordinates.
  That is, the complex (real and imaginary) will have radius (magnitude) r
  and angle (phase) phi.

  #### See also
  [new/2](#new/2) [imag/0](#imag/0)

  #### Examples
      iex> Complex.from_polar(1, :math.pi/2)
      %Complex{im: 1.0, re: 6.123233995736766e-17}

      iex> polar = Complex.to_polar(%Complex{re: 6, im: 20})
      iex> Complex.from_polar(polar)
      %Complex{im: 20.0, re: 6.0}
  """
  @spec from_polar({number, number}) :: complex
  def from_polar({r, phi}), do: from_polar(r, phi)

  @spec from_polar(number, number) :: complex
  def from_polar(r, phi) do
    new(r * :math.cos(phi), r * :math.sin(phi))
  end

  @doc """
  Returns the phase angle of the supplied complex, in radians.

  #### See also
  [new/2](#new/2) [from_polar/2](#from_polar/2)

  #### Examples
      iex> Complex.phase( Complex.from_polar(1,:math.pi/2) )
      1.5707963267948966
  """
  @spec phase(complex) :: float
  @spec phase(number) :: number
  def phase(n) when n < 0, do: :math.pi()
  def phase(n) when is_number(n), do: 0

  def phase(z = %Complex{}) do
    :math.atan2(z.im, z.re)
  end

  @doc """
  Returns the polar coordinates of the supplied complex.  That is, the
  returned tuple {r,phi} is the magnitude and phase (in radians) of z.

  #### See also
  [from_polar/2](#from_polar/2)

  #### Examples
      iex> Complex.to_polar( Complex.from_polar(1,:math.pi/2) )
      {1.0, 1.5707963267948966}
  """
  @spec to_polar(complex) :: {float, float}
  def to_polar(complex)

  def to_polar(z = %Complex{}) do
    {abs(z), phase(z)}
  end

  def to_polar(n) when n < 0, do: {n, :math.pi()}

  def to_polar(n), do: {n, 0}

  @doc """
    Returns a new complex that is the sum of the provided complex numbers.  Also
    supports a mix of complex and number.

    #### See also
    [div/2](#div/2), [multiply/2](#multiply/2), [subtract/2](#subtract/2)

    #### Examples
        iex> Complex.add( Complex.from_polar(1, :math.pi/2), Complex.from_polar(1, :math.pi/2) )
        %Complex{im: 2.0, re: 1.2246467991473532e-16}

        iex> Complex.add( Complex.new(4, 4), 1 )
        %Complex{im: 4, re: 5}

        iex> Complex.add( 2, Complex.new(4, 3) )
        %Complex{im: 3, re: 6}

        iex> Complex.add(2, 3)
        5

        iex> Complex.add(2.0, 2)
        4.0
  """
  @spec add(complex, complex) :: complex
  @spec add(number, complex) :: complex
  @spec add(complex, number) :: complex
  @spec add(number, number) :: number
  def add(left, right) when is_number(left) and is_number(right), do: left + right

  def add(left, right) do
    %Complex{re: re_left, im: im_left} = as_complex(left)
    %Complex{re: re_right, im: im_right} = as_complex(right)
    new(re_left + re_right, im_left + im_right)
  end

  @doc """
    Returns a new complex that is the difference of the provided complex numbers.
    Also supports a mix of complex and number.

    #### See also
    [add/2](#add/2), [div/2](#div/2), [multiply/2](#multiply/2)

    #### Examples
        iex> Complex.subtract( Complex.new(1,2), Complex.new(3,4) )
        %Complex{im: -2, re: 0-2}

        iex> Complex.subtract( Complex.new(1, 2), 3 )
        %Complex{im: 2, re: -2}

        iex> Complex.subtract( 10, Complex.new(1, 2) )
        %Complex{im: -2, re: 9}
  """
  @spec subtract(complex, complex) :: complex
  @spec subtract(number, complex) :: complex
  @spec subtract(complex, number) :: complex
  @spec subtract(number, number) :: number
  def subtract(left, right) when is_number(left) and is_number(right), do: left - right

  def subtract(left, right) do
    %Complex{re: re_left, im: im_left} = as_complex(left)
    %Complex{re: re_right, im: im_right} = as_complex(right)
    new(re_left - re_right, im_left - im_right)
  end

  @doc """
  Returns a new complex that is the product of the provided complex numbers.
  Also supports a mix of complex and number.

  #### See also
  [add/2](#add/2), [div/2](#div/2), [subtract/2](#subtract/2)

  #### Examples
      iex> Complex.multiply( Complex.new(1,2), Complex.new(3,4) )
      %Complex{im: 10, re: -5}

      iex> Complex.multiply( Complex.imag(), Complex.imag() )
      %Complex{im: 0.0, re: -1.0}

      iex> Complex.multiply(Complex.new(1, 2), 3 )
      %Complex{im: 6, re: 3}

      iex> Complex.multiply( 3, Complex.new(1, 2) )
      %Complex{im: 6, re: 3}
  """
  @spec multiply(complex, complex) :: complex
  @spec multiply(number, complex) :: complex
  @spec multiply(complex, number) :: complex
  @spec multiply(number, number) :: number
  def multiply(left, right) when is_number(left) and is_number(right), do: left * right

  def multiply(left, right) do
    %Complex{re: r1, im: i1} = as_complex(left)
    %Complex{re: r2, im: i2} = as_complex(right)
    new(r1 * r2 - i1 * i2, i1 * r2 + r1 * i2)
  end

  @doc """
  Returns a new complex that is the square of the provided complex number.

  #### See also
  [multiply/2](#multiply/2)

  #### Examples
      iex> Complex.square( Complex.new(2.0, 0.0) )
      %Complex{im: 0.0, re: 4.0}

      iex> Complex.square( Complex.imag() )
      %Complex{im: 0.0, re: -1.0}
  """
  @spec square(complex) :: complex
  @spec square(number) :: number
  def square(z), do: multiply(z, z)

  @doc """
  Returns a new complex that is the ratio (division) of the provided complex
  numbers.

  #### See also
  [add/2](#add/2), [multiply/2](#multiply/2), [subtract/2](#subtract/2)

  #### Examples
      iex> Complex.divide( Complex.from_polar(1, :math.pi/2), Complex.from_polar(1, :math.pi/2) )
      %Complex{im: 0.0, re: 1.0}
  """
  @spec divide(complex, complex) :: complex
  @spec divide(number, complex) :: complex
  @spec divide(complex, number) :: complex
  @spec divide(number, number) :: number

  def divide(x, y) when is_number(x) and is_number(y), do: x / y

  def divide(x, y) do
    %Complex{re: r1, im: i1} = as_complex(x)
    %Complex{re: r2, im: i2} = as_complex(y)

    if Kernel.abs(r2) < Kernel.abs(i2) do
      r = r2 / i2
      den = i2 + r * r2
      new((r1 * r + i1) / den, (i1 * r - r1) / den)
    else
      r = i2 / r2
      den = r2 + r * i2
      new((r1 + r * i1) / den, (i1 - r * r1) / den)
    end
  end

  @doc """
  Returns the magnitude (length) of the provided complex number.

  #### See also
  [new/2](#new/2), [phase/1](#phase/1)

  #### Examples
      iex> Complex.abs( Complex.from_polar(1, :math.pi/2) )
      1.0
  """
  @spec abs(complex) :: number
  @spec abs(number) :: number
  def abs(n) when is_number(n), do: Kernel.abs(n)

  def abs(%Complex{re: r, im: i}) do
    # optimized by checking special cases (sqrt is expensive)
    x = Kernel.abs(r)
    y = Kernel.abs(i)

    cond do
      x == 0.0 -> y
      y == 0.0 -> x
      x > y -> x * :math.sqrt(1.0 + y / x * (y / x))
      true -> y * :math.sqrt(1.0 + x / y * (x / y))
    end
  end

  @doc """
  Returns the square of the magnitude of the provided complex number.

  The square of the magnitude is faster to compute---no square roots!

  #### See also
  [new/2](#new/2), [abs/1](#abs/1)

  #### Examples
      iex> Complex.abs_squared( Complex.from_polar(1, :math.pi/2) )
      1.0

      iex> Complex.abs_squared( Complex.from_polar(2, :math.pi/2) )
      4.0
  """
  @spec abs_squared(complex) :: number
  @spec abs_squared(number) :: number
  def abs_squared(n) when is_number(n), do: n * n

  def abs_squared(%Complex{re: r, im: i}) do
    r * r + i * i
  end

  @doc """
  Returns a new complex that is the complex conjugate of the provided complex
  number.

  #### See also
  [abs/2](#abs/2), [phase/1](#phase/1)

  #### Examples
      iex> Complex.conjugate( Complex.new(1,2) )
      %Complex{im: -2, re: 1}
  """
  @spec conjugate(complex) :: complex
  @spec conjugate(number) :: number
  def conjugate(n) when is_number(n), do: n

  def conjugate(%Complex{re: r, im: i}) do
    new(r, -i)
  end

  @doc """
  Returns a new complex that is the complex square root of the provided
  complex number.

  #### See also
  [abs/2](#abs/2), [phase/1](#phase/1)

  #### Examples
      iex> Complex.sqrt( Complex.from_polar(2,:math.pi) )
      %Complex{im: 1.4142135623730951, re: 8.659560562354933e-17}
  """
  @spec sqrt(complex) :: complex
  @spec sqrt(number) :: number
  def sqrt(n) when is_number(n), do: :math.sqrt(n)

  def sqrt(z = %Complex{re: r, im: i}) do
    if z.re == 0.0 and z.im == 0.0 do
      new(z.re, z.im)
    else
      x = Kernel.abs(r)
      y = Kernel.abs(i)

      w =
        if x >= y do
          :math.sqrt(x) * :math.sqrt(0.5 * (1.0 + :math.sqrt(1.0 + y / x * (y / x))))
        else
          :math.sqrt(y) * :math.sqrt(0.5 * (x / y + :math.sqrt(1.0 + x / y * (x / y))))
        end

      if z.re >= 0.0 do
        new(w, z.im / (2 * w))
      else
        i2 =
          if z.im >= 0.0 do
            w
          else
            -w
          end

        new(z.im / (2 * i2), i2)
      end
    end
  end

  @doc """
  Returns a new complex that is the complex exponential of the provided
  complex number.  That is, e raised to the power z.

  #### See also
  [ln/1](#ln/1)

  #### Examples
      iex> Complex.exp( Complex.from_polar(2,:math.pi) )
      %Complex{im: 3.3147584285483636e-17, re: 0.1353352832366127}
  """
  @spec exp(complex) :: complex
  @spec exp(number) :: number
  def exp(n) when is_number(n), do: :math.exp(n)

  def exp(z = %Complex{}) do
    rho = :math.exp(z.re)
    theta = z.im
    new(rho * :math.cos(theta), rho * :math.sin(theta))
  end

  @doc """
  Returns a new complex that is the complex natural log of the provided
  complex number.  That is, log base e of z.

  #### See also
  [exp/1](#exp/1)

  #### Examples
      iex> Complex.ln( Complex.from_polar(2,:math.pi) )
      %Complex{im: 3.141592653589793, re: 0.6931471805599453}
  """
  @spec ln(complex) :: complex
  def ln(n) when is_number(n), do: :math.log(n)

  def ln(z = %Complex{}) do
    new(:math.log(abs(z)), :math.atan2(z.im, z.re))
  end

  @doc """
  Returns a new complex that is the complex log base 10 of the provided
  complex number.

  #### See also
  [ln/1](#ln/1)

  #### Examples
      iex> Complex.log10( Complex.from_polar(2,:math.pi) )
      %Complex{im: 1.3643763538418412, re: 0.30102999566398114}
  """
  @spec log10(complex) :: complex
  @spec log10(number) :: number

  def log10(n) when is_number(n), do: :math.log10(n)

  def log10(z = %Complex{}) do
    divide(ln(z), new(:math.log(10.0), 0.0))
  end

  @doc """
  Returns a new complex that is the complex log base 2 of the provided
  complex number.

  #### See also
  [ln/1](#ln/1), [log10/1](#log10/1)

  #### Examples
      iex> Complex.log2( Complex.from_polar(2,:math.pi) )
      %Complex{im: 4.532360141827194, re: 1.0}
  """
  @spec log2(complex) :: complex
  @spec log2(number) :: number

  def log2(n) when is_number(n), do: :math.log2(n)

  def log2(z = %Complex{}) do
    divide(ln(z), new(:math.log(2.0), 0.0))
  end

  @doc """
  Returns a new complex that is the provided parameter a raised to the
  complex power b.

  #### See also
  [ln/1](#ln/1), [log10/1](#log10/1)

  #### Examples
      iex> Complex.power( Complex.from_polar(2,:math.pi), Complex.imag() )
      %Complex{im: 0.027612020368333014, re: 0.03324182700885666}
  """
  @spec power(complex, complex) :: complex
  @spec power(number, complex) :: complex
  @spec power(complex, number) :: complex
  @spec power(number, number) :: number

  def power(x, y) when is_number(x) and is_number(y), do: :math.pow(x, y)

  def power(x, y) do
    x = as_complex(x)
    y = as_complex(y)

    cond do
      x.re == 0.0 and x.im == 0.0 ->
        if y.re == 0.0 and y.im == 0.0 do
          new(1.0, 0.0)
        else
          new(0.0, 0.0)
        end

      y.re == 1.0 and y.im == 0.0 ->
        x

      y.re == -1.0 and y.im == 0.0 ->
        divide(new(1.0, 0.0), x)

      true ->
        rho = :math.sqrt(x.re * x.re + x.im * x.im)
        theta = :math.atan2(x.im, x.re)
        s = :math.pow(rho, y.re) * :math.exp(-y.im * theta)
        r = y.re * theta + y.im * :math.log(rho)
        new(s * :math.cos(r), s * :math.sin(r))
    end
  end

  @doc """
  Returns a new complex that is the sine of the provided parameter.

  #### See also
  [cos/1](#cos/1), [tan/1](#tan/1)

  #### Examples
      iex> Complex.sin( Complex.from_polar(2,:math.pi) )
      %Complex{im: -1.0192657827055095e-16, re: -0.9092974268256817}
  """
  @spec sin(complex) :: complex
  @spec sin(number) :: number

  def sin(n) when is_number(n), do: :math.sin(n)

  def sin(z = %Complex{}) do
    new(
      :math.sin(z.re) * :math.cosh(z.im),
      :math.cos(z.re) * :math.sinh(z.im)
    )
  end

  @doc """
  Returns a new complex that is the "negation" of the provided parameter.
  That is, the real and imaginary parts are negated.

  #### See also
  [neq/2](#new/2), [imag/0](#imag/0)

  #### Examples
      iex> Complex.negate( Complex.new(3,5) )
      %Complex{im: -5, re: -3}
  """
  @spec negate(complex) :: complex
  @spec negate(number) :: number
  def negate(n) when is_number(n), do: -n

  def negate(z = %Complex{}) do
    new(-z.re, -z.im)
  end

  @doc """
  Returns a new complex that is the inverse sine (i.e., arcsine) of the
  provided parameter.

  #### See also
  [sin/1](#sin/1)

  #### Examples
      iex> Complex.asin( Complex.from_polar(2,:math.pi) )
      %Complex{im: 1.3169578969248164, re: -1.5707963267948966}
  """
  @spec asin(complex) :: complex
  @spec asin(number) :: number
  def asin(n) when is_number(n), do: :math.asin(n)

  def asin(z = %Complex{}) do
    i = new(0.0, 1.0)
    # result = -i*ln(i*z + sqrt(1.0-z*z))
    # result = -i*ln(t1 + sqrt(t2))
    t1 = multiply(i, z)
    t2 = subtract(new(1.0, 0.0), multiply(z, z))
    multiply(negate(i), ln(add(t1, sqrt(t2))))
  end

  @doc """
  Returns a new complex that is the cosine of the provided parameter.

  #### See also
  [sin/1](#sin/1), [tan/1](#tan/1)

  #### Examples
      iex> Complex.cos( Complex.from_polar(2,:math.pi) )
      %Complex{im: 2.2271363664699914e-16, re: -0.4161468365471424}
  """
  @spec cos(complex) :: complex
  @spec cos(number) :: number
  def cos(n) when is_number(n), do: :math.cos(n)

  def cos(z = %Complex{}) do
    new(
      :math.cos(z.re) * :math.cosh(z.im),
      -:math.sin(z.re) * :math.sinh(z.im)
    )
  end

  @doc """
  Returns a new complex that is the inverse cosine (i.e., arccosine) of the
  provided parameter.

  #### See also
  [cos/1](#cos/1)

  #### Examples
      iex> Complex.acos( Complex.from_polar(2,:math.pi) )
      %Complex{im: 1.3169578969248164, re: -3.141592653589793}
  """
  @spec acos(complex) :: complex
  @spec acos(number) :: number
  def acos(n) when is_number(n), do: :math.acos(n)

  def acos(z = %Complex{}) do
    i = new(0.0, 1.0)
    one = new(1.0, 0.0)
    # result = -i*ln(z + sqrt(z*z-1.0))
    # result = -i*ln(z + sqrt(t1))
    t1 = subtract(multiply(z, z), one)
    multiply(negate(i), ln(add(z, sqrt(t1))))
  end

  @doc """
  Returns a new complex that is the tangent of the provided parameter.

  #### See also
  [sin/1](#sin/1), [cos/1](#cos/1)

  #### Examples
      iex> Complex.tan( Complex.from_polar(2,:math.pi) )
      %Complex{im: 1.4143199004457917e-15, re: 2.185039863261519}
  """
  @spec tan(complex) :: complex
  @spec tan(number) :: number
  def tan(n) when is_number(n), do: :math.tan(n)

  def tan(z = %Complex{}) do
    divide(sin(z), cos(z))
  end

  @doc """
  Returns a new complex that is the inverse tangent (i.e., arctangent) of the
  provided parameter.

  #### See also
  [tan/1](#tan/1), [atan2/2](#atan2/2)

  #### Examples
      iex> Complex.atan( Complex.from_polar(2,:math.pi) )
      %Complex{im: 0.0, re: -1.1071487177940904}

      iex> Complex.tan( Complex.atan(Complex.new(2,3)) )
      %Complex{im: 3.0, re: 2.0}
  """
  @spec atan(complex) :: complex
  @spec atan(number) :: number
  def atan(n) when is_number(n), do: :math.atan(n)

  def atan(z = %Complex{}) do
    i = new(0.0, 1.0)
    # result = 0.5*i*(ln(1-i*z)-ln(1+i*z))
    t1 = multiply(new(0.5, 0.0), i)
    t2 = subtract(new(1.0, 0.0), multiply(i, z))
    t3 = add(new(1.0, 0.0), multiply(i, z))
    multiply(t1, subtract(ln(t2), ln(t3)))
  end

  @doc """
  Returns the phase of the complex number `b + i*a`.

  #### See also
  [tan/1](#tan/1), [atan/1](#atan/1)

  #### Examples
      iex> phase = Complex.atan2(2, 2)
      iex> phase == :math.pi() / 4
      true

      iex> phase = Complex.atan2(2, Complex.new(0))
      iex> phase == Complex.new(:math.pi() / 2, 0)
      true
  """
  def atan2(a, b) when is_number(a) and is_number(b), do: :math.atan2(a, b)

  def atan2(a, b) do
    a = as_complex(a)
    b = as_complex(b)

    if a.im != 0 or b.im != 0 do
      raise ArithmeticError, "Complex.atan2 only accepts real numbers as arguments"
    end

    a.re
    |> :math.atan2(b.re)
    |> Complex.new()
  end

  @doc """
  Returns a new complex that is the cotangent of the provided parameter.

  #### See also
  [sin/1](#sin/1), [cos/1](#cos/1), [tan/1](#tan/1)

  #### Examples
      iex> Complex.cot( Complex.from_polar(2,:math.pi) )
      %Complex{im: -2.9622992129532336e-16, re: 0.45765755436028577}
  """
  @spec cot(complex) :: complex
  @spec cot(number) :: number
  def cot(n) when is_number(n), do: 1 / :math.tan(n)

  def cot(z = %Complex{}) do
    divide(cos(z), sin(z))
  end

  @doc """
  Returns a new complex that is the inverse cotangent (i.e., arccotangent) of
  the provided parameter.

  #### See also
  [cot/1](#cot/1)

  #### Examples
      iex> Complex.acot( Complex.from_polar(2,:math.pi) )
      %Complex{im: -9.71445146547012e-17, re: -0.46364760900080615}

      iex> Complex.cot( Complex.acot(Complex.new(2,3)) )
      %Complex{im: 3.0, re: 1.9999999999999991}
  """
  @spec acot(complex) :: complex
  @spec acot(number) :: number
  def acot(n) when is_number(n), do: :math.atan(1 / n)

  def acot(z = %Complex{}) do
    i = new(0.0, 1.0)
    # result = 0.5*i*(ln(1-i/z)-ln(1+i/z))
    t1 = multiply(new(0.5, 0.0), i)
    t2 = subtract(new(1.0, 0.0), divide(i, z))
    t3 = add(new(1.0, 0.0), divide(i, z))
    multiply(t1, subtract(ln(t2), ln(t3)))
  end

  @doc """
  Returns a new complex that is the secant of the provided parameter.

  #### See also
  [sin/1](#sin/1), [cos/1](#cos/1), [tan/1](#tan/1)

  #### Examples
      iex> Complex.sec( Complex.from_polar(2,:math.pi) )
      %Complex{im: -1.2860374461837126e-15, re: -2.402997961722381}
  """
  @spec sec(complex) :: complex
  @spec sec(number) :: number
  def sec(z) do
    divide(1, cos(z))
  end

  @doc """
  Returns a new complex that is the inverse secant (i.e., arcsecant) of
  the provided parameter.

  #### See also
  [sec/1](#sec/1)

  #### Examples
      iex> Complex.asec( Complex.from_polar(2,:math.pi) )
      %Complex{im: 0.0, re: 2.0943951023931957}

      iex> Complex.sec( Complex.asec(Complex.new(2,3)) )
      %Complex{im: 2.9999999999999982, re: 1.9999999999999987}
  """
  @spec asec(complex) :: complex
  @spec asec(number) :: number
  def asec(n) when is_number(n) do
    :math.acos(1 / n)
  end

  def asec(z = %Complex{}) do
    i = new(0.0, 1.0)
    # result = -i*ln(i*sqrt(1-1/(z*z))+1/z)
    # result = -i*ln(i*sqrt(1-t2)+t1)
    t1 = divide(1, z)
    t2 = square(t1)
    # result = -i*ln(i*sqrt(t3)+t1)
    # result = -i*ln(t4+t1)
    t3 = subtract(1, t2)
    t4 = multiply(i, sqrt(t3))
    multiply(negate(i), ln(add(t4, t1)))
  end

  @doc """
  Returns a new complex that is the cosecant of the provided parameter.

  #### See also
  [sec/1](#sec/1), [sin/1](#sin/1), [cos/1](#cos/1), [tan/1](#tan/1)

  #### Examples
      iex> Complex.csc( Complex.from_polar(2,:math.pi) )
      %Complex{im: 1.2327514463765779e-16, re: -1.0997501702946164}
  """
  @spec csc(complex) :: complex
  @spec csc(number) :: number
  def csc(z) do
    divide(1, sin(z))
  end

  @doc """
  Returns a new complex that is the inverse cosecant (i.e., arccosecant) of
  the provided parameter.

  #### See also
  [sec/1](#sec/1)

  #### Examples
      iex> Complex.acsc( Complex.from_polar(2,:math.pi) )
      %Complex{im: 0.0, re: -0.5235987755982988}

      iex> Complex.csc( Complex.acsc(Complex.new(2,3)) )
      %Complex{im: 2.9999999999999996, re: 1.9999999999999991}
  """
  @spec acsc(complex) :: complex
  @spec acsc(number) :: number
  def acsc(n) when is_number(n), do: :math.asin(1 / n)

  def acsc(z = %Complex{}) do
    i = new(0.0, 1.0)
    one = new(1.0, 0.0)
    # result = -i*ln(sqrt(1-1/(z*z))+i/z)
    # result = -i*ln(sqrt(1-t2)+t1)
    t1 = divide(i, z)
    t2 = divide(one, multiply(z, z))
    # result = -i*ln(sqrt(t3)+t1)
    # result = -i*ln(t4+t1)
    t3 = subtract(one, t2)
    t4 = sqrt(t3)
    multiply(negate(i), ln(add(t4, t1)))
  end

  @doc """
  Returns a new complex that is the hyperbolic sine of the provided parameter.

  #### See also
  [cosh/1](#cosh/1), [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.sinh( Complex.from_polar(2,:math.pi) )
      %Complex{im: 9.214721821703068e-16, re: -3.626860407847019}
  """
  @spec sinh(complex) :: complex
  @spec sinh(number) :: number
  def sinh(n) when is_number(n), do: :math.sinh(n)

  def sinh(z = %Complex{}) do
    z
    |> exp()
    |> subtract(exp(negate(z)))
    |> divide(2)
  end

  @doc """
  Returns a new complex that is the inverse hyperbolic sine (i.e., arcsinh) of
  the provided parameter.

  #### See also
  [sinh/1](#sinh/1)

  #### Examples
      iex> Complex.asinh( Complex.from_polar(2,:math.pi) )
      %Complex{im: 1.0953573965284052e-16, re: -1.4436354751788099}

      iex> Complex.sinh( Complex.asinh(Complex.new(2,3)) )
      %Complex{im: 3.0, re: 2.0000000000000004}
  """
  @spec asinh(complex) :: complex
  @spec asinh(number) :: number

  if math_fun_supported?(:asinh, 1) do
    def asinh(n) when is_number(n) do
      :math.asinh(n)
    end
  else
    def asinh(n) when is_number(n) do
      :math.log(n + :math.sqrt(1 + n * n))
    end
  end

  def asinh(z) do
    # result = ln(z+sqrt(z*z+1))
    # result = ln(z+sqrt(t1))
    # result = ln(t2)
    t1 = add(multiply(z, z), 1)
    t2 = add(z, sqrt(t1))
    ln(t2)
  end

  @doc """
  Returns a new complex that is the hyperbolic cosine of the provided
  parameter.

  #### See also
  [sinh/1](#sinh/1), [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.cosh( Complex.from_polar(2,:math.pi) )
      %Complex{im: -8.883245978848233e-16, re: 3.7621956910836314}
  """
  @spec cosh(complex) :: complex
  @spec cosh(number) :: number
  def cosh(z) do
    z
    |> exp()
    |> add(exp(negate(z)))
    |> divide(2)
  end

  @doc """
  Returns a new complex that is the inverse hyperbolic cosine (i.e., arccosh)
  of the provided parameter.

  #### See also
  [cosh/1](#cosh/1)

  #### Examples
      iex> Complex.acosh( Complex.from_polar(2,:math.pi) )
      %Complex{im: -3.141592653589793, re: -1.3169578969248164}
  """
  # Windows apparently doesn't support :math.acosh
  @spec acosh(complex) :: complex
  @spec acosh(number) :: number
  if math_fun_supported?(:acosh, 1) do
    def acosh(n) when is_number(n), do: :math.acosh(n)
  else
    def acosh(n) when is_number(n) do
      :math.log(n + :math.sqrt(n * n - 1))
    end
  end

  def acosh(z = %Complex{}) do
    # result = ln(z+sqrt(z*z-1))
    # result = ln(z+sqrt(t1))
    # result = ln(t2)
    t1 = subtract(multiply(z, z), 1)
    t2 = add(z, sqrt(t1))
    ln(t2)
  end

  @doc """
  Returns a new complex that is the hyperbolic tangent of the provided
  parameter.

  #### See also
  [sinh/1](#sinh/1), [cosh/1](#cosh/1)

  #### Examples
      iex> Complex.tanh( Complex.from_polar(2,:math.pi) )
      %Complex{im: 1.7304461302709572e-17, re: -0.964027580075817}
  """
  @spec tanh(complex) :: complex
  @spec tanh(number) :: number
  def tanh(n) when is_number(n), do: :math.tanh(n)

  def tanh(z = %Complex{}) do
    divide(sinh(z), cosh(z))
  end

  @doc """
  Returns a new complex that is the inverse hyperbolic tangent (i.e., arctanh)
  of the provided parameter.

  #### See also
  [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.atanh( Complex.from_polar(2,:math.pi) )
      %Complex{im: 1.5707963267948966, re: -0.5493061443340549}

      iex> Complex.tanh( Complex.atanh(Complex.new(2,3)) )
      %Complex{im: 2.999999999999999, re: 1.9999999999999987}
  """
  @spec atanh(complex) :: complex
  @spec atanh(number) :: number
  if math_fun_supported?(:atanh, 1) do
    def atanh(n) when is_number(n), do: :math.atanh(n)
  else
    def atanh(n) when is_number(n) do
      0.5 * :math.log((1 + n) / (1 - n))
    end
  end

  def atanh(z = %Complex{}) do
    one = new(1.0, 0.0)
    p5 = new(0.5, 0.0)
    # result = 0.5*(ln((1+z)/(1-z)))
    # result = 0.5*(ln(t2/t1))
    # result = 0.5*(ln(t3))
    t1 = subtract(one, z)
    t2 = add(one, z)
    t3 = divide(t2, t1)
    multiply(p5, ln(t3))
  end

  @doc """
  Returns a new complex that is the hyperbolic secant of the provided
  parameter.

  #### See also
  [sinh/1](#sinh/1), [cosh/1](#cosh/1), [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.sech( Complex.from_polar(2,:math.pi) )
      %Complex{im: 6.27608655779184e-17, re: 0.2658022288340797}
  """
  @spec sech(complex) :: complex
  @spec sech(number) :: number
  def sech(z) do
    divide(1, cosh(z))
  end

  @doc """
  Returns a new complex that is the inverse hyperbolic secant (i.e., arcsech)
  of the provided parameter.

  #### See also
  [sech/1](#sech/1)

  #### Examples
      iex> Complex.asech( Complex.from_polar(2,:math.pi) )
      %Complex{im: -2.0943951023931953, re: 0.0}

      iex> Complex.sech( Complex.asech(Complex.new(2,3)) )
      %Complex{im: 2.999999999999999, re: 2.0}
  """
  @spec asech(complex) :: complex
  @spec asech(number) :: number
  def asech(z) do
    # result = ln(1/z+sqrt(1/z+1)*sqrt(1/z-1))
    # result = ln(t1+sqrt(t1+1)*sqrt(t1-1))
    # result = ln(t1+t2*t3)
    t1 = divide(1, z)
    t2 = sqrt(add(t1, 1))
    t3 = sqrt(subtract(t1, 1))
    ln(add(t1, multiply(t2, t3)))
  end

  @doc """
  Returns a new complex that is the hyperbolic cosecant of the provided
  parameter.

  #### See also
  [sinh/1](#sinh/1), [cosh/1](#cosh/1), [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.csch( Complex.from_polar(2,:math.pi) )
      %Complex{im: -7.00520014334671e-17, re: -0.2757205647717832}
  """
  @spec csch(complex) :: complex
  @spec csch(number) :: number
  def csch(z), do: divide(1, sinh(z))

  @doc """
  Returns a new complex that is the inverse hyperbolic cosecant (i.e., arccsch)
  of the provided parameter.

  #### See also
  [csch/1](#csch/1)

  #### Examples
      iex> Complex.acsch( Complex.from_polar(2,:math.pi) )
      %Complex{im: -5.4767869826420256e-17, re: -0.48121182505960336}

      iex> Complex.csch( Complex.acsch(Complex.new(2,3)) )
      %Complex{im: 3.0000000000000018, re: 1.9999999999999982}
  """
  @spec acsch(complex) :: complex
  @spec acsch(number) :: number
  def acsch(z) do
    # result = ln(1/z+sqrt(1/(z*z)+1))
    # result = ln(t1+sqrt(t2+1))
    # result = ln(t1+t3)
    t1 = divide(1, z)
    t2 = divide(1, multiply(z, z))
    t3 = sqrt(add(t2, 1))
    ln(add(t1, t3))
  end

  @doc """
  Returns a new complex that is the hyperbolic cotangent of the provided
  parameter.

  #### See also
  [sinh/1](#sinh/1), [cosh/1](#cosh/1), [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.coth( Complex.from_polar(2,:math.pi) )
      %Complex{im: -1.8619978115303632e-17, re: -1.037314720727548}
  """
  @spec coth(complex) :: complex
  @spec coth(number) :: number
  def coth(z) do
    divide(cosh(z), sinh(z))
  end

  @doc """
  Returns a new complex that is the inverse hyperbolic cotangent (i.e., arccoth)
  of the provided parameter.

  #### See also
  [coth/1](#coth/1)

  #### Examples
      iex> Complex.acoth( Complex.from_polar(2,:math.pi) )
      %Complex{im: -8.164311994315688e-17, re: -0.5493061443340548}

      iex> Complex.coth( Complex.acoth(Complex.new(2,3)) )
      %Complex{im: 2.999999999999998, re: 2.000000000000001}
  """
  @spec acoth(complex) :: complex
  @spec acoth(number) :: number
  def acoth(z) do
    # result = 0.5*(ln(1+1/z)-ln(1-1/z))
    # result = 0.5*(ln(1+t1)-ln(1-t1))
    # result = 0.5*(ln(t2)-ln(t3))
    t1 = divide(1, z)
    t2 = add(1, t1)
    t3 = subtract(1, t1)
    multiply(0.5, subtract(ln(t2), ln(t3)))
  end

  @doc """
  Calculates erf(x) of the argument
  """

  if math_fun_supported?(:erf, 1) do
    def erf(x) when is_number(x) do
      :math.erf(x)
    end
  else
    def erf(x) when is_number(x) do
      x = x |> max(-4.0) |> min(4.0)
      x2 = x * x

      alpha =
        0.0
        |> muladd(x2, -2.72614225801306e-10)
        |> muladd(x2, 2.77068142495902e-08)
        |> muladd(x2, -2.10102402082508e-06)
        |> muladd(x2, -5.69250639462346e-05)
        |> muladd(x2, -7.34990630326855e-04)
        |> muladd(x2, -2.95459980854025e-03)
        |> muladd(x2, -1.60960333262415e-02)

      beta =
        0.0
        |> muladd(x2, -1.45660718464996e-05)
        |> muladd(x2, -2.13374055278905e-04)
        |> muladd(x2, -1.68282697438203e-03)
        |> muladd(x2, -7.37332916720468e-03)
        |> muladd(x2, -1.42647390514189e-02)

      min(x * alpha / beta, 1.0)
    end
  end

  def erf(%Complex{re: re, im: im}) do
    if im != 0 do
      raise ArithmeticError, "erf not implemented for non-real numbers"
    end

    Complex.new(erf(re))
  end

  if math_fun_supported?(:erfc, 1) do
    def erfc(z) when is_number(z), do: :math.erfc(z)
  end

  def erfc(z), do: subtract(1, erf(z))

  def erf_inv(z) when is_number(z) do
    w = -:math.log((1 - z) * (1 + z))
    erf_inv_p(w) * z
  end

  def erf_inv(%Complex{re: re, im: im}) do
    if im != 0 do
      raise ArithmeticError, "erf_inv not implemented for non-real numbers"
    end

    Complex.new(erf_inv(re))
  end

  defp erf_inv_p(w) when w < 5 do
    w = w - 2.5

    2.81022636e-08
    |> muladd(w, 3.43273939e-07)
    |> muladd(w, -3.5233877e-06)
    |> muladd(w, -4.39150654e-06)
    |> muladd(w, 0.00021858087)
    |> muladd(w, -0.00125372503)
    |> muladd(w, -0.00417768164)
    |> muladd(w, 0.246640727)
    |> muladd(w, 1.50140941)
  end

  defp erf_inv_p(w) do
    w = :math.sqrt(w) - 3

    -0.000200214257
    |> muladd(w, 0.000100950558)
    |> muladd(w, 0.00134934322)
    |> muladd(w, -0.00367342844)
    |> muladd(w, 0.00573950773)
    |> muladd(w, -0.0076224613)
    |> muladd(w, 0.00943887047)
    |> muladd(w, 1.00167406)
    |> muladd(w, 2.83297682)
  end

  defp muladd(acc, t, n) do
    acc * t + n
  end

  ### Logical Operations
  @doc "Logical equality between two numbers"
  def equal?(a, b), do: a == b

  @doc "Logical difference between two numbers"
  def not_equal?(a, b), do: a != b

  @doc """
  Logical AND between two numbers.

  Zero is false, everything else is true.
  """
  def logical_and?(a, b), do: as_boolean(a) and as_boolean(b)

  @doc """
  Logical OR between two numbers.

  Zero is false, everything else is true.
  """
  def logical_or?(a, b), do: as_boolean(a) or as_boolean(b)

  @doc """
  Logical XOR between two numbers.

  Zero is false, everything else is true.
  """
  def logical_xor?(a, b), do: as_boolean(a) != as_boolean(b)

  @doc """
  Component-wise `>` between two complex numbers.
  """
  def greater?(a, b) when is_number(a) and is_number(b), do: a > b

  def greater?(a, b) do
    %{re: a_re, im: a_im} = as_complex(a)
    %{re: b_re, im: b_im} = as_complex(b)

    a_re > b_re and a_im > b_im
  end

  @doc """
  Component-wise `<` between two complex numbers.
  """
  def less?(a, b) when is_number(a) and is_number(b), do: a < b

  def less?(a, b) do
    %{re: a_re, im: a_im} = as_complex(a)
    %{re: b_re, im: b_im} = as_complex(b)

    a_re < b_re and a_im < b_im
  end

  @doc """
  Component-wise `>=` between two complex numbers.
  """
  def greater_equal?(a, b) when is_number(a) and is_number(b), do: a >= b

  def greater_equal?(a, b) do
    %{re: a_re, im: a_im} = as_complex(a)
    %{re: b_re, im: b_im} = as_complex(b)

    a_re >= b_re and a_im >= b_im
  end

  @doc """
  Component-wise `<=` between two complex numbers.
  """
  def less_equal?(a, b) when is_number(a) and is_number(b), do: a <= b

  def less_equal?(a, b) do
    %{re: a_re, im: a_im} = as_complex(a)
    %{re: b_re, im: b_im} = as_complex(b)

    a_re <= b_re and a_im <= b_im
  end

  defp as_complex(%Complex{} = x), do: x
  defp as_complex(x) when is_number(x), do: new(x)

  defp as_boolean(n) when n == 0, do: false
  defp as_boolean(%Complex{re: re, im: im}) when re == 0 and im == 0, do: false
  defp as_boolean(_), do: true
end
