defmodule Complex do
  @moduledoc """
  *Complex* is a library for types and mathematical functions for complex
  numbers.

  Each complex number is represented as a structure holding the real and
  imaginary part.  There are functions for creation and manipulation of
  them.  Unfortunately since there is no operator overloading in Elixir the
  math functions (add, subtract, etc.) are implemented as add/2, sub/2, etc.

  ## Examples
      iex> Complex.new(3, 4)
      %Complex{im: 4, re: 3}

      iex> Complex.imag()
      %Complex{im: 1.0, re: 0.0}
  """

  @vsn 2

  import Kernel, except: [abs: 1, div: 2]

  @typedoc """
  General type for complex numbers
  """
  @type complex :: %Complex{re: number, im: number}
  defstruct re: 0, im: 0

  # Distinguishes between a Complex type and a number type.  This function
  # is used to allow mixed Complex and numbers in function calls.
  defp decompose(number) when is_number(number), do: {number, 0}
  defp decompose(%Complex{re: re, im: im}), do: {re, im}

  @doc """
  Returns a new complex with specified real and imaginary components.  The
  imaginary part defaults to zero so a "real" number can be created with new/1

  #### See also
  [imag/0](#imag/0) [fromPolar/2](#fromPolar/2)

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
  [new/2](#new/2) [fromPolar/2](#fromPolar/2)

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
      iex> Complex.fromPolar(1, :math.pi/2)
      %Complex{im: 1.0, re: 6.123233995736766e-17}
  """
  @spec fromPolar(number, number) :: complex
  def fromPolar(r, phi) do
    new(r * :math.cos(phi), r * :math.sin(phi))
  end

  @doc """
  Returns the phase angle of the supplied complex, in radians.

  #### See also
  [new/2](#new/2) [fromPolar/2](#fromPolar/2)

  #### Examples
      iex> Complex.phase( Complex.fromPolar(1,:math.pi/2) )
      1.5707963267948966
  """
  @spec phase(complex) :: float
  def phase(z = %Complex{}) do
    :math.atan2(z.im, z.re)
  end

  @doc """
  Returns the polar coordinates of the supplied complex.  That is, the
  returned tuple {r,phi} is the magnitude and phase (in radians) of z.

  #### See also
  [fromPolar/2](#fromPolar/2)

  #### Examples
      iex> Complex.getPolar( Complex.fromPolar(1,:math.pi/2) )
      {1.0, 1.5707963267948966}
  """
  @spec getPolar(complex) :: {float, float}
  def getPolar(z = %Complex{}) do
    {abs(z), phase(z)}
  end

  @doc """
    Returns a new complex that is the sum of the provided complex numbers.  Also
    supports a mix of complex and number.

    #### See also
    [div/2](#div/2), [mult/2](#mult/2), [sub/2](#sub/2)

    #### Examples
        iex> Complex.add( Complex.fromPolar(1, :math.pi/2), Complex.fromPolar(1, :math.pi/2) )
        %Complex{im: 2.0, re: 1.2246467991473532e-16}

        iex> Complex.add( Complex.new(4, 4), 1 )
        %Complex{im: 4, re: 5}

        iex> Complex.add( 2, Complex.new(4, 3) )
        %Complex{im: 3, re: 6}

        iex> Complex.add( 2, 3 )
        %Complex{im: 0, re: 5}
  """
  @spec add(complex, complex) :: complex
  @spec add(number, complex) :: complex
  @spec add(complex, number) :: complex
  def add(left, right) when is_number(left) and is_number(right), do: new(left + right, 0)

  def add(left, right) do
    {re_left, im_left} = decompose(left)
    {re_right, im_right} = decompose(right)
    new(re_left + re_right, im_left + im_right)
  end

  @doc """
    Returns a new complex that is the difference of the provided complex numbers.
    Also supports a mix of complex and number.

    #### See also
    [add/2](#add/2), [div/2](#div/2), [mult/2](#mult/2)

    #### Examples
        iex> Complex.sub( Complex.new(1,2), Complex.new(3,4) )
        %Complex{im: -2, re: 0-2}

        iex> Complex.sub( Complex.new(1, 2), 3 )
        %Complex{im: 2, re: -2}

        iex> Complex.sub( 10, Complex.new(1, 2) )
        %Complex{im: -2, re: 9}
  """
  @spec sub(complex, complex) :: complex
  @spec sub(number, complex) :: complex
  @spec sub(complex, number) :: complex
  def sub(left, right) when is_number(left) and is_number(right), do: new(left - right, 0)

  def sub(left, right) do
    {re_left, im_left} = decompose(left)
    {re_right, im_right} = decompose(right)
    new(re_left - re_right, im_left - im_right)
  end

  @doc """
  Returns a new complex that is the product of the provided complex numbers.
  Also supports a mix of complex and number.

  #### See also
  [add/2](#add/2), [div/2](#div/2), [sub/2](#sub/2)

  #### Examples
      iex> Complex.mult( Complex.new(1,2), Complex.new(3,4) )
      %Complex{im: 10, re: -5}

      iex> Complex.mult( Complex.imag(), Complex.imag() )
      %Complex{im: 0.0, re: -1.0}

      iex> Complex.mult(Complex.new(1, 2), 3 )
      %Complex{im: 6, re: 3}

      iex> Complex.mult( 3, Complex.new(1, 2) )
      %Complex{im: 6, re: 3}
  """
  @spec mult(complex, complex) :: complex
  @spec mult(number, complex) :: complex
  @spec mult(complex, number) :: complex
  def mult(left, right) when is_number(left) and is_number(right), do: new(left * right, 0)

  def mult(left, right) do
    {r1, i1} = decompose(left)
    {r2, i2} = decompose(right)
    new(r1 * r2 - i1 * i2, i1 * r2 + r1 * i2)
  end

  @doc """
  Returns a new complex that is the square of the provided complex number.

  #### See also
  [mult/2](#mult/2)

  #### Examples
      iex> Complex.square( Complex.new(2.0, 0.0) )
      %Complex{im: 0.0, re: 4.0}

      iex> Complex.square( Complex.imag() )
      %Complex{im: 0.0, re: -1.0}
  """
  @spec square(complex) :: complex
  def square(z), do: mult(z, z)

  @doc """
  Returns a new complex that is the ratio (division) of the provided complex
  numbers.

  #### See also
  [add/2](#add/2), [mult/2](#mult/2), [sub/2](#sub/2)

  #### Examples
      iex> Complex.div( Complex.fromPolar(1, :math.pi/2), Complex.fromPolar(1, :math.pi/2) )
      %Complex{im: 0.0, re: 1.0}
  """
  @spec div(complex, complex) :: complex
  def div(%Complex{re: r1, im: i1}, %Complex{re: r2, im: i2}) do
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
      iex> Complex.abs( Complex.fromPolar(1, :math.pi/2) )
      1.0
  """
  @spec abs(complex) :: number
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
      iex> Complex.abs_squared( Complex.fromPolar(1, :math.pi/2) )
      1.0

      iex> Complex.abs_squared( Complex.fromPolar(2, :math.pi/2) )
      4.0
  """
  @spec abs_squared(complex) :: number
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
  def conjugate(%Complex{re: r, im: i}) do
    new(r, -i)
  end

  @doc """
  Returns a new complex that is the complex square root of the provided
  complex number.

  #### See also
  [abs/2](#abs/2), [phase/1](#phase/1)

  #### Examples
      iex> Complex.sqrt( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 1.4142135623730951, re: 8.659560562354933e-17}
  """
  @spec sqrt(complex) :: complex
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
      iex> Complex.exp( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 3.3147584285483636e-17, re: 0.1353352832366127}
  """
  @spec exp(complex) :: complex
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
      iex> Complex.ln( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 3.141592653589793, re: 0.6931471805599453}
  """
  @spec ln(complex) :: complex
  def ln(z = %Complex{}) do
    new(:math.log(abs(z)), :math.atan2(z.im, z.re))
  end

  @doc """
  Returns a new complex that is the complex log base 10 of the provided
  complex number.

  #### See also
  [ln/1](#ln/1)

  #### Examples
      iex> Complex.log10( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 1.3643763538418412, re: 0.30102999566398114}
  """
  @spec log10(complex) :: complex
  def log10(z = %Complex{}) do
    div(ln(z), new(:math.log(10.0), 0.0))
  end

  @doc """
  Returns a new complex that is the complex log base 2 of the provided
  complex number.

  #### See also
  [ln/1](#ln/1), [log10/1](#log10/1)

  #### Examples
      iex> Complex.log2( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 4.532360141827194, re: 1.0}
  """
  @spec log2(complex) :: complex
  def log2(z = %Complex{}) do
    div(ln(z), new(:math.log(2.0), 0.0))
  end

  @doc """
  Returns a new complex that is the provided parameter a raised to the
  complex power b.

  #### See also
  [ln/1](#ln/1), [log10/1](#log10/1)

  #### Examples
      iex> Complex.pow( Complex.fromPolar(2,:math.pi), Complex.imag() )
      %Complex{im: 0.027612020368333014, re: 0.03324182700885666}
  """
  @spec pow(complex, complex) :: complex
  def pow(x = %Complex{}, y = %Complex{}) do
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
        div(new(1.0, 0.0), x)

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
      iex> Complex.sin( Complex.fromPolar(2,:math.pi) )
      %Complex{im: -1.0192657827055095e-16, re: -0.9092974268256817}
  """
  @spec sin(complex) :: complex
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
      iex> Complex.neg( Complex.new(3,5) )
      %Complex{im: -5, re: -3}
  """
  @spec neg(complex) :: complex
  def neg(z = %Complex{}) do
    new(-z.re, -z.im)
  end

  @doc """
  Returns a new complex that is the inverse sine (i.e., arcsine) of the
  provided parameter.

  #### See also
  [sin/1](#sin/1)

  #### Examples
      iex> Complex.asin( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 1.3169578969248164, re: -1.5707963267948963}

      iex> Complex.sin( Complex.asin(Complex.new(2,3)) )
      %Complex{im: 3.000000000000001, re: 1.9999999999999991}
  """
  @spec asin(complex) :: complex
  def asin(z = %Complex{}) do
    i = new(0.0, 1.0)
    # result = -i*ln(i*z + sqrt(1.0-z*z))
    # result = -i*ln(t1 + sqrt(t2))
    t1 = mult(i, z)
    t2 = sub(new(1.0, 0.0), mult(z, z))
    mult(neg(i), ln(add(t1, sqrt(t2))))
  end

  @doc """
  Returns a new complex that is the cosine of the provided parameter.

  #### See also
  [sin/1](#sin/1), [tan/1](#tan/1)

  #### Examples
      iex> Complex.cos( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 2.2271363664699914e-16, re: -0.4161468365471424}
  """
  @spec cos(complex) :: complex
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
      iex> Complex.acos( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 1.3169578969248164, re: -3.141592653589793}

      iex> Complex.cos( Complex.acos(Complex.new(2,3)) )
      %Complex{im: 3.0, re: 2.0000000000000004}
  """
  @spec acos(complex) :: complex
  def acos(z = %Complex{}) do
    i = new(0.0, 1.0)
    one = new(1.0, 0.0)
    # result = -i*ln(z + sqrt(z*z-1.0))
    # result = -i*ln(z + sqrt(t1))
    t1 = sub(mult(z, z), one)
    mult(neg(i), ln(add(z, sqrt(t1))))
  end

  @doc """
  Returns a new complex that is the tangent of the provided parameter.

  #### See also
  [sin/1](#sin/1), [cos/1](#cos/1)

  #### Examples
      iex> Complex.tan( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 1.4143199004457917e-15, re: 2.185039863261519}
  """
  @spec tan(complex) :: complex
  def tan(z = %Complex{}) do
    div(sin(z), cos(z))
  end

  @doc """
  Returns a new complex that is the inverse tangent (i.e., arctangent) of the
  provided parameter.

  #### See also
  [tan/1](#tan/1)

  #### Examples
      iex> Complex.atan( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 0.0, re: -1.1071487177940904}

      iex> Complex.tan( Complex.atan(Complex.new(2,3)) )
      %Complex{im: 3.0, re: 2.0}
  """
  @spec atan(complex) :: complex
  def atan(z = %Complex{}) do
    i = new(0.0, 1.0)
    # result = 0.5*i*(ln(1-i*z)-ln(1+i*z))
    t1 = mult(new(0.5, 0.0), i)
    t2 = sub(new(1.0, 0.0), mult(i, z))
    t3 = add(new(1.0, 0.0), mult(i, z))
    mult(t1, sub(ln(t2), ln(t3)))
  end

  @doc """
  Returns a new complex that is the cotangent of the provided parameter.

  #### See also
  [sin/1](#sin/1), [cos/1](#cos/1), [tan/1](#tan/1)

  #### Examples
      iex> Complex.cot( Complex.fromPolar(2,:math.pi) )
      %Complex{im: -2.9622992129532336e-16, re: 0.45765755436028577}
  """
  @spec cot(complex) :: complex
  def cot(z = %Complex{}) do
    div(cos(z), sin(z))
  end

  @doc """
  Returns a new complex that is the inverse cotangent (i.e., arccotangent) of
  the provided parameter.

  #### See also
  [cot/1](#cot/1)

  #### Examples
      iex> Complex.acot( Complex.fromPolar(2,:math.pi) )
      %Complex{im: -9.71445146547012e-17, re: -0.46364760900080615}

      iex> Complex.cot( Complex.acot(Complex.new(2,3)) )
      %Complex{im: 2.9999999999999996, re: 1.9999999999999991}
  """
  @spec acot(complex) :: complex
  def acot(z = %Complex{}) do
    i = new(0.0, 1.0)
    # result = 0.5*i*(ln(1-i/z)-ln(1+i/z))
    t1 = mult(new(0.5, 0.0), i)
    t2 = sub(new(1.0, 0.0), div(i, z))
    t3 = add(new(1.0, 0.0), div(i, z))
    mult(t1, sub(ln(t2), ln(t3)))
  end

  @doc """
  Returns a new complex that is the secant of the provided parameter.

  #### See also
  [sin/1](#sin/1), [cos/1](#cos/1), [tan/1](#tan/1)

  #### Examples
      iex> Complex.sec( Complex.fromPolar(2,:math.pi) )
      %Complex{im: -1.2860374461837126e-15, re: -2.402997961722381}
  """
  @spec sec(complex) :: complex
  def sec(z = %Complex{}) do
    div(new(1.0, 0.0), cos(z))
  end

  @doc """
  Returns a new complex that is the inverse secant (i.e., arcsecant) of
  the provided parameter.

  #### See also
  [sec/1](#sec/1)

  #### Examples
      iex> Complex.asec( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 0.0, re: 2.0943951023931957}

      iex> Complex.sec( Complex.asec(Complex.new(2,3)) )
      %Complex{im: 2.9999999999999987, re: 1.9999999999999987}
  """
  @spec asec(complex) :: complex
  def asec(z = %Complex{}) do
    i = new(0.0, 1.0)
    one = new(1.0, 0.0)
    # result = -i*ln(i*sqrt(1-1/(z*z))+1/z)
    # result = -i*ln(i*sqrt(1-t2)+t1)
    t1 = div(one, z)
    t2 = div(one, mult(z, z))
    # result = -i*ln(i*sqrt(t3)+t1)
    # result = -i*ln(t4+t1)
    t3 = sub(one, t2)
    t4 = mult(i, sqrt(t3))
    mult(neg(i), ln(add(t4, t1)))
  end

  @doc """
  Returns a new complex that is the cosecant of the provided parameter.

  #### See also
  [sec/1](#sec/1), [sin/1](#sin/1), [cos/1](#cos/1), [tan/1](#tan/1)

  #### Examples
      iex> Complex.csc( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 1.2327514463765779e-16, re: -1.0997501702946164}
  """
  @spec csc(complex) :: complex
  def csc(z = %Complex{}) do
    one = new(1.0, 0.0)
    div(one, sin(z))
  end

  @doc """
  Returns a new complex that is the inverse cosecant (i.e., arccosecant) of
  the provided parameter.

  #### See also
  [sec/1](#sec/1)

  #### Examples
      iex> Complex.acsc( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 0.0, re: -0.5235987755982988}

      iex> Complex.csc( Complex.acsc(Complex.new(2,3)) )
      %Complex{im: 2.9999999999999996, re: 1.9999999999999993}
  """
  @spec acsc(complex) :: complex
  def acsc(z = %Complex{}) do
    i = new(0.0, 1.0)
    one = new(1.0, 0.0)
    # result = -i*ln(sqrt(1-1/(z*z))+i/z)
    # result = -i*ln(sqrt(1-t2)+t1)
    t1 = div(i, z)
    t2 = div(one, mult(z, z))
    # result = -i*ln(sqrt(t3)+t1)
    # result = -i*ln(t4+t1)
    t3 = sub(one, t2)
    t4 = sqrt(t3)
    mult(neg(i), ln(add(t4, t1)))
  end

  @doc """
  Returns a new complex that is the hyperbolic sine of the provided parameter.

  #### See also
  [cosh/1](#cosh/1), [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.sinh( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 9.214721821703068e-16, re: -3.626860407847019}
  """
  @spec sinh(complex) :: complex
  def sinh(z = %Complex{}) do
    p5 = new(0.5, 0.0)
    mult(p5, sub(exp(z), exp(neg(z))))
  end

  @doc """
  Returns a new complex that is the inverse hyperbolic sine (i.e., arcsinh) of
  the provided parameter.

  #### See also
  [sinh/1](#sinh/1)

  #### Examples
      iex> Complex.asinh( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 1.0953573965284052e-16, re: -1.4436354751788099}

      iex> Complex.sinh( Complex.asinh(Complex.new(2,3)) )
      %Complex{im: 3.0, re: 2.000000000000001}
  """
  @spec asinh(complex) :: complex
  def asinh(z = %Complex{}) do
    one = new(1.0, 0.0)
    # result = ln(z+sqrt(z*z+1))
    # result = ln(z+sqrt(t1))
    # result = ln(t2)
    t1 = add(mult(z, z), one)
    t2 = add(z, sqrt(t1))
    ln(t2)
  end

  @doc """
  Returns a new complex that is the hyperbolic cosine of the provided
  parameter.

  #### See also
  [sinh/1](#sinh/1), [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.cosh( Complex.fromPolar(2,:math.pi) )
      %Complex{im: -8.883245978848233e-16, re: 3.7621956910836314}
  """
  @spec cosh(complex) :: complex
  def cosh(z = %Complex{}) do
    p5 = new(0.5, 0.0)
    mult(p5, add(exp(z), exp(neg(z))))
  end

  @doc """
  Returns a new complex that is the inverse hyperbolic cosine (i.e., arccosh)
  of the provided parameter.

  #### See also
  [cosh/1](#cosh/1)

  #### Examples
      iex> Complex.acosh( Complex.fromPolar(2,:math.pi) )
      %Complex{im: -3.141592653589793, re: -1.3169578969248164}

      iex> Complex.cosh( Complex.acosh(Complex.new(2,3)) )
      %Complex{im: 3.0, re: 2.0}
  """
  @spec acosh(complex) :: complex
  def acosh(z = %Complex{}) do
    one = new(1.0, 0.0)
    # result = ln(z+sqrt(z*z-1))
    # result = ln(z+sqrt(t1))
    # result = ln(t2)
    t1 = sub(mult(z, z), one)
    t2 = add(z, sqrt(t1))
    ln(t2)
  end

  @doc """
  Returns a new complex that is the hyperbolic tangent of the provided
  parameter.

  #### See also
  [sinh/1](#sinh/1), [cosh/1](#cosh/1)

  #### Examples
      iex> Complex.tanh( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 1.7304461302709572e-17, re: -0.964027580075817}
  """
  @spec tanh(complex) :: complex
  def tanh(z = %Complex{}) do
    div(sinh(z), cosh(z))
  end

  @doc """
  Returns a new complex that is the inverse hyperbolic tangent (i.e., arctanh)
  of the provided parameter.

  #### See also
  [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.atanh( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 1.5707963267948966, re: -0.5493061443340549}

      iex> Complex.tanh( Complex.atanh(Complex.new(2,3)) )
      %Complex{im: 2.999999999999999, re: 1.9999999999999987}
  """
  @spec atanh(complex) :: complex
  def atanh(z = %Complex{}) do
    one = new(1.0, 0.0)
    p5 = new(0.5, 0.0)
    # result = 0.5*(ln((1+z)/(1-z)))
    # result = 0.5*(ln(t2/t1))
    # result = 0.5*(ln(t3))
    t1 = sub(one, z)
    t2 = add(one, z)
    t3 = div(t2, t1)
    mult(p5, ln(t3))
  end

  @doc """
  Returns a new complex that is the hyperbolic secant of the provided
  parameter.

  #### See also
  [sinh/1](#sinh/1), [cosh/1](#cosh/1), [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.sech( Complex.fromPolar(2,:math.pi) )
      %Complex{im: 6.27608655779184e-17, re: 0.2658022288340797}
  """
  @spec sech(complex) :: complex
  def sech(z = %Complex{}) do
    two = new(2.0, 0.0)
    div(two, add(exp(z), exp(neg(z))))
  end

  @doc """
  Returns a new complex that is the inverse hyperbolic secant (i.e., arcsech)
  of the provided parameter.

  #### See also
  [sech/1](#sech/1)

  #### Examples
      iex> Complex.asech( Complex.fromPolar(2,:math.pi) )
      %Complex{im: -2.0943951023931953, re: 0.0}

      iex> Complex.sech( Complex.asech(Complex.new(2,3)) )
      %Complex{im: 2.999999999999999, re: 2.0}
  """
  @spec asech(complex) :: complex
  def asech(z = %Complex{}) do
    one = new(1.0, 0.0)
    # result = ln(1/z+sqrt(1/z+1)*sqrt(1/z-1))
    # result = ln(t1+sqrt(t1+1)*sqrt(t1-1))
    # result = ln(t1+t2*t3)
    t1 = div(one, z)
    t2 = sqrt(add(t1, one))
    t3 = sqrt(sub(t1, one))
    ln(add(t1, mult(t2, t3)))
  end

  @doc """
  Returns a new complex that is the hyperbolic cosecant of the provided
  parameter.

  #### See also
  [sinh/1](#sinh/1), [cosh/1](#cosh/1), [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.csch( Complex.fromPolar(2,:math.pi) )
      %Complex{im: -7.00520014334671e-17, re: -0.2757205647717832}
  """
  @spec csch(complex) :: complex
  def csch(z = %Complex{}) do
    two = new(2.0, 0.0)
    div(two, sub(exp(z), exp(neg(z))))
  end

  @doc """
  Returns a new complex that is the inverse hyperbolic cosecant (i.e., arccsch)
  of the provided parameter.

  #### See also
  [csch/1](#csch/1)

  #### Examples
      iex> Complex.acsch( Complex.fromPolar(2,:math.pi) )
      %Complex{im: -5.4767869826420256e-17, re: -0.48121182505960336}

      iex> Complex.csch( Complex.acsch(Complex.new(2,3)) )
      %Complex{im: 3.0000000000000018, re: 1.9999999999999982}
  """
  @spec acsch(complex) :: complex
  def acsch(z = %Complex{}) do
    one = new(1.0, 0.0)
    # result = ln(1/z+sqrt(1/(z*z)+1))
    # result = ln(t1+sqrt(t2+1))
    # result = ln(t1+t3)
    t1 = div(one, z)
    t2 = div(one, mult(z, z))
    t3 = sqrt(add(t2, one))
    ln(add(t1, t3))
  end

  @doc """
  Returns a new complex that is the hyperbolic cotangent of the provided
  parameter.

  #### See also
  [sinh/1](#sinh/1), [cosh/1](#cosh/1), [tanh/1](#tanh/1)

  #### Examples
      iex> Complex.coth( Complex.fromPolar(2,:math.pi) )
      %Complex{im: -1.8619978115303632e-17, re: -1.037314720727548}
  """
  @spec coth(complex) :: complex
  def coth(z = %Complex{}) do
    div(cosh(z), sinh(z))
  end

  @doc """
  Returns a new complex that is the inverse hyperbolic cotangent (i.e., arccoth)
  of the provided parameter.

  #### See also
  [coth/1](#coth/1)

  #### Examples
      iex> Complex.acoth( Complex.fromPolar(2,:math.pi) )
      %Complex{im: -8.164311994315688e-17, re: -0.5493061443340548}

      iex> Complex.coth( Complex.acoth(Complex.new(2,3)) )
      %Complex{im: 2.999999999999998, re: 2.000000000000001}
  """
  @spec acoth(complex) :: complex
  def acoth(z = %Complex{}) do
    one = new(1.0, 0.0)
    p5 = new(0.5, 0.0)
    # result = 0.5*(ln(1+1/z)-ln(1-1/z))
    # result = 0.5*(ln(1+t1)-ln(1-t1))
    # result = 0.5*(ln(t2)-ln(t3))
    t1 = div(one, z)
    t2 = add(one, t1)
    t3 = sub(one, t1)
    mult(p5, sub(ln(t2), ln(t3)))
  end
end
