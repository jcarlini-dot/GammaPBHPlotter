import sys, argparse
import numpy as np
import matplotlib.pyplot as plt
import gammapbh as gp

def main():
    p = argparse.ArgumentParser(description="Plot PBH gamma-ray spectra for a single mass.")
    p.add_argument("--mass", type=float, required=True, help="BH mass (must match a grid entry)")
    p.add_argument("--align", choices=["secondary","primary","union","intersection"], default="secondary")
    args = p.parse_args()

    E, comps = gp.load_spectra_components(args.mass, align=args.align)
    y = comps.get("total")
    if y is None:
        sys.exit("No 'total' spectrum exposed by API.")
    plt.figure()
    plt.loglog(E, y, label=f"total ({args.align})")
    plt.xlabel("Energy")
    plt.ylabel("dN/dE (arb.)")
    plt.title(f"PBH spectra, M = {args.mass:.3e}")
    plt.legend()
    plt.grid(True, which="both", ls=":")
    plt.show()

if __name__ == "__main__":
    main()