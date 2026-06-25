import sys
import argparse
from haplo_bench.generator import generate

def main():
    parser = argparse.ArgumentParser(description="HAPLO-BENCH CLI tool")
    subparsers = parser.add_subparsers(dest="command")
    
    # generate command
    gen_parser = subparsers.add_parser("generate", help="Generate benchmark instance")
    gen_parser.add_argument("--config", required=True, help="Path to config YAML")
    gen_parser.add_argument("--out", required=True, help="Output directory")
    gen_parser.add_argument("--seed", type=int, help="Random seed")
    
    args = parser.parse_args()
    
    if args.command == "generate":
        generate(args.config, args.out, args.seed)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
