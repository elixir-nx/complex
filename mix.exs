defmodule Complex.Mixfile do
  use Mix.Project

  def project do
    [
      app: :complex,
      version: "0.3.0",
      description: description(),
      package: package(),
      elixir: "~> 1.1",
      deps: deps(),
      docs: [extras: []]
    ]
  end

  # Configuration for the OTP application
  #
  # Type "mix help compile.app" for more information
  def application do
    [applications: [:logger]]
  end

  defp deps do
    [
      {:ex_doc, "~> 0.22", only: :dev, runtime: false},
      {:dialyxir, "~> 0.4", only: [:dev]}
    ]
  end

  defp description do
    """
    Complex is a library for types and mathematical functions for complex
    numbers.
    """
  end

  defp package do
    [
      maintainers: ["Tom Krauss"],
      licenses: ["Apache 2.0"],
      links: %{
        "GitHub" => "https://github.com/twist-vector/elixir-complex.git",
        "Docs" => "http://hexdocs.pm/complex"
      }
    ]
  end
end
