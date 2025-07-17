#!/usr/bin/env python3
"""
Generate Fortran macros for MetStateType fields.

This script parses a Fortran source file containing the MetStateType definition and generates Fortran include files
for field accessors, allocation, and deallocation macros. The generated macros are used to implement field-by-field
allocation and access in a modular, robust way.

Usage
-----
Run from the command line:
    generate_metstate_macros.py accessor <input_f90> <output_inc>
    generate_metstate_macros.py allocate <input_f90> <output_inc>
    generate_metstate_macros.py deallocate <input_f90> <output_inc>
    generate_metstate_macros.py column_accessor <input_f90> <output_inc>
    generate_metstate_macros.py 2d_scalar_accessor <input_f90> <output_inc>
    generate_metstate_macros.py scalar_accessor <input_f90> <output_inc>

Parameters
----------
mode : str
    Which macro to generate. One of: accessor, allocate, deallocate, column_accessor, 2d_scalar_accessor, scalar_accessor.
input_f90 : str
    Path to the Fortran source file containing the MetStateType definition.
output_inc : str
    Path to the output .inc file to write.

Examples
--------
$ python generate_metstate_macros.py allocate ../src/core/metstate_mod.F90 metstate_allocate_fields.inc

"""
import re
import sys

def parse_metstate_type(filename):
    """
    Parse the MetStateType definition in a Fortran source file.

    Parameters
    ----------
    filename : str
        Path to the Fortran source file containing MetStateType.

    Returns
    -------
    fields : list of tuple
        List of (name, rank, dims) for each allocatable REAL(fp) field.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    in_type = False
    fields = []  # List of (name, rank, dims)
    for line in lines:
        if 'TYPE, PUBLIC :: MetStateType' in line:
            in_type = True
        elif in_type and 'end type' in line.lower():
            break
        elif in_type:
            # Match real(fp), allocatable :: name(dimensions)
            m = re.match(r'\s*REAL\(fp\),\s*ALLOCATABLE\s*::\s*(\w+)\s*\(([^)]*)\)', line, re.IGNORECASE)
            if m:
                name = m.group(1)
                dims = m.group(2)
                rank = dims.count(',') + 1
                fields.append((name, rank, dims))

    return fields

def write_accessor(fields, output_file):
    """
    Write a macro for field accessors (all 2D/3D fields).

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    # Deduplicate by case-insensitive field name
    unique_fields = {}
    for name, rank, dims in fields:
        key = name.lower()
        if key not in unique_fields:
            unique_fields[key] = (name, rank, dims)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, rank, dims = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case(" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) then\n")
            if rank == 3:
                f.write(f"      column_ptr => this%{name}(col_i, col_j, :)\n")
            elif rank == 2:
                f.write(f"      column_ptr => null()\n")
            f.write(f"   endif\n")
        f.write("case default\n")
        f.write("   column_ptr => null()\n")

def write_allocate(fields, output_file):
    """
    Write a macro for allocation statements for MetStateType fields.
    Skips scalar (rank-0) fields. Adds a case for 'ALL'.

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    with open(output_file, 'w') as f:
        # ALL case
        f.write("case ('ALL', 'all')\n")
        for name, rank, dims in fields:
            if rank == 0:
                continue
            if name.lower() in ("soilm"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoil))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            elif name.lower() in ("frsoil"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoiltype))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            elif name.lower() in ("frlanduse", "frlai"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nSURFTYPE))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            else:
                if rank == 2:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
                elif rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nz))\n")
                else:
                    raise ValueError(f"Unsupported rank {rank} for field {name}")
        # Individual field cases
        for name, rank, dims in fields:
            if rank == 0:
                continue  # Skip scalars
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            if name.lower() in ("soilm"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoil))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            elif name.lower() in ("frsoil"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nsoiltype))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            elif name.lower() in ("frlanduse", "frlai"):
                if rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nSURFTYPE))\n")
                else:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
            else:
                if rank == 2:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny))\n")
                elif rank == 3:
                    f.write(f"  if (.not.allocated(this%{name})) allocate(this%{name}(nx,ny,nz))\n")
                else:
                    raise ValueError(f"Unsupported rank {rank} for field {name}")
        f.write("case default\n")
        f.write("  ! No allocation for unknown field\n")

def write_deallocate(fields, output_file):
    """
    Write a macro for deallocation statements for MetStateType fields.
    Skips scalar (rank-0) fields. Adds a case for 'ALL'.

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    with open(output_file, 'w') as f:
        # ALL case
        f.write("case ('ALL', 'all')\n")
        for name, rank, dims in fields:
            if rank == 0:
                continue
            f.write(f"  if (allocated(this%{name})) deallocate(this%{name})\n")
        # Individual field cases
        for name, rank, dims in fields:
            if rank == 0:
                continue  # Skip scalars
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"  if (allocated(this%{name})) deallocate(this%{name})\n")
        f.write("case default\n")
        f.write("  ! No deallocation for unknown field\n")

def write_column_accessor(fields, output_file):
    """
    Write a macro for 3D column field accessors.

    Parameters
    ----------
    fields : list of tuple
        List of (name, rank, dims) for each field.
    output_file : str
        Path to output .inc file.
    """
    # Only 3D fields
    unique_fields = {}
    for name, rank, dims in fields:
        if rank == 3:
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, rank, dims)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, rank, dims = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) column_ptr => this%{name}(col_i, col_j, :)\n")
        f.write("case default\n")
        f.write("   column_ptr => null()\n")

def write_2d_scalar_accessor(fields, output_file):
    """
    Write a macro for 2D scalar field accessors.

    Parameters
    ----------
    fields : list of (name, rank, dims)
        List of fields.
    output_file : str
        Path to output .inc file.
    """
    # Only 2D fields
    unique_fields = {}
    for name, rank, dims in fields:
        if rank == 2:
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, rank, dims)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, rank, dims = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   if (allocated(this%{name})) scalar_val = this%{name}(col_i, col_j)\n")
        f.write("case default\n")
        f.write("   scalar_val = 0.0_fp\n")

def write_scalar_accessor(fields, output_file):
    """
    Write a macro for scalar (rank-0) field accessors.

    Parameters
    ----------
    fields : list of (name, rank, dims)
        List of fields.
    output_file : str
        Path to output .inc file.
    """
    # Only scalar fields (rank 0)
    unique_fields = {}
    for name, rank, dims in fields:
        if rank == 0:
            key = name.lower()
            if key not in unique_fields:
                unique_fields[key] = (name, rank, dims)
    with open(output_file, 'w') as f:
        for key in sorted(unique_fields):
            name, rank, dims = unique_fields[key]
            labels = sorted({name, name.lower()})
            f.write("case (" + ", ".join(f"'{label}'" for label in labels) + ")\n")
            f.write(f"   scalar_val = this%{name}\n")
        f.write("case default\n")
        f.write("   scalar_val = 0.0_fp\n")

def main():
    """
    Main entry point for the macro generator script.

    Parses command-line arguments and dispatches to the appropriate macro writer.
    """
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate Fortran macros for MetStateType fields.",
        epilog="Example: generate_metstate_macros.py allocate metstate_mod.F90 metstate_allocate_fields.inc"
    )
    parser.add_argument('mode', choices=['accessor', 'allocate', 'deallocate', 'column_accessor', '2d_scalar_accessor', 'scalar_accessor'],
                        help="Type of macro to generate.")
    parser.add_argument('input_f90', help="Input Fortran file containing MetStateType.")
    parser.add_argument('output_inc', help="Output .inc file to write.")
    args = parser.parse_args()

    fields = parse_metstate_type(args.input_f90)
    if not fields:
        print(f"Warning: No fields found in {args.input_f90}.")

    if args.mode == 'accessor':
        write_accessor(fields, args.output_inc)
    elif args.mode == 'allocate':
        write_allocate(fields, args.output_inc)
    elif args.mode == 'deallocate':
        write_deallocate(fields, args.output_inc)
    elif args.mode == 'column_accessor':
        write_column_accessor(fields, args.output_inc)
    elif args.mode == '2d_scalar_accessor':
        write_2d_scalar_accessor(fields, args.output_inc)
    elif args.mode == 'scalar_accessor':
        write_scalar_accessor(fields, args.output_inc)
    else:
        parser.error(f"Unknown mode: {args.mode}")

if __name__ == "__main__":
    main()
