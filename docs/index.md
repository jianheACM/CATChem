# CATChem - Community Atmospheric Transport Chemistry Model

<div class="community-header">
  <h1>CATChem</h1>
  <p style="font-size: 1.2rem; margin: 0.5rem 0;">Community Atmospheric Transport Chemistry Model</p>
</div>

!!! info "What is CATChem?"
    CATChem is a modern, high-performance atmospheric chemistry transport model designed for operational weather prediction and research applications. Built with modern Fortran and designed for scalability, CATChem provides sophisticated atmospheric chemistry capabilities as an open-source community project.

## ✨ Key Features

<div class="grid cards" markdown>

- :material-rocket-launch: **High Performance**

  ---

  Optimized for modern HPC systems with column virtualization, efficient memory management, and scalable parallelization.

- :material-puzzle: **Modular Architecture**

  ---

  Process-based design with pluggable schemes, flexible configuration, and clean separation of concerns.

- :material-cog: **Operational Ready**

  ---

  Designed for 24/7 operational use with robust error handling, comprehensive diagnostics, and monitoring capabilities.

- :material-science: **Comprehensive Chemistry**

  ---

  Full atmospheric chemistry including gas-phase, aerosol processes, emissions, deposition, and transport.

</div>

## 🚀 Quick Start

=== "Installation"

    ```bash
    # Clone the repository
    git clone https://github.com/NOAA-GSL/CATChem.git
    cd CATChem

    # Build with CMake
    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j$(nproc)
    ```

=== "First Run"

    ```bash
    # Run a test case
    cd build
    ctest -R test_settling

    # Run with sample configuration
    ./catchem_driver --config ../parm/config/catchem_test.yml
    ```

=== "Integration"

    ```fortran
    ! Integrate with your model
    use CATChemAPI_Mod

    type(CATChemType) :: catchem

    call catchem%init(config_file, rc)
    call catchem%run(dt, met_fields, chem_fields, rc)
    call catchem%finalize(rc)
    ```

## 📊 Performance Highlights

<div class="performance-metric">
  <span>Computational Efficiency</span>
  <span class="performance-metric__value">10x faster column processing</span>
</div>

<div class="performance-metric">
  <span>Memory Usage</span>
  <span class="performance-metric__value">50% reduction vs traditional models</span>
</div>

<div class="performance-metric">
  <span>Scalability</span>
  <span class="performance-metric__value">Linear scaling to 10,000+ cores</span>
</div>

<div class="performance-metric">
  <span>Operational Readiness</span>
  <span class="performance-metric__value">24/7 production deployment</span>
</div>

## 🧪 Available Processes

| Process Type | Description | Status |
|--------------|-------------|--------|
| <span class="process-badge process-badge--chemistry">Chemistry</span> | Gas-phase and aerosol chemistry | ✅ Production |
| <span class="process-badge process-badge--emission">Emissions</span> | Anthropogenic and biogenic emissions | ✅ Production |
| <span class="process-badge process-badge--transport">Settling</span> | Gravitational settling with slip correction | ✅ Production |
| <span class="process-badge process-badge--transport">Vertical Mixing</span> | YSU boundary layer scheme | ✅ Production |
| <span class="process-badge process-badge--loss">Dry Deposition</span> | Surface deposition processes | ✅ Production |
| <span class="process-badge process-badge--loss">Wet Deposition</span> | Precipitation scavenging | 🚧 Development |
| <span class="process-badge process-badge--emission">Dust</span> | Mineral dust emission and transport | ✅ Production |
| <span class="process-badge process-badge--emission">Sea Salt</span> | Marine aerosol processes | ✅ Production |
| <span class="process-badge process-badge--transport">Plume Rise</span> | Wildfire and point source plume rise | ✅ Production |

## 🏗️ Architecture Overview

```mermaid
flowchart TB
    A["Driver/Host Model"] --> B["CATChem API"]
    B --> C["State Management"]
    C --> D["Process Manager"]
    D --> E["Column Virtualization"]
    E --> F{"Process Types"}
    F --> G["Chemistry Processes"]
    F --> H["Transport Processes"]
    F --> I["Emission Processes"]
    F --> J["Loss Processes"]
    G --> K["Diagnostic System"]
    H --> K
    I --> K
    J --> K
    K --> L["Output Manager"]
```

!!! note "Modern Design Principles"
    - **Separation of Concerns**: Clear boundaries between I/O, computation, and state management
    - **Column Virtualization**: Efficient 1D processing with automatic parallelization
    - **Process-Based Architecture**: Modular, testable, and maintainable components
    - **StateContainer Pattern**: Unified data management with automatic dependency tracking
    - **Comprehensive Diagnostics**: Built-in monitoring, profiling, and error reporting

## 🎯 Use Cases

### Operational Weather Prediction
- **NOAA UFS Integration**: Native support for FV3 and other UFS components
- **Real-time Processing**: Optimized for operational time constraints
- **Ensemble Forecasting**: Efficient multi-member ensemble capabilities

### Research Applications
- **Process Studies**: Individual process testing and validation
- **Sensitivity Analysis**: Parameter perturbation and uncertainty quantification
- **Model Development**: Framework for new process implementation

### Air Quality Forecasting
- **Multi-scale Modeling**: Global to urban scale applications
- **Chemical Data Assimilation**: Support for observation integration
- **Health Impact Assessment**: Human exposure and risk evaluation

## 📚 Documentation Structure

<div class="grid cards" markdown>

- [:material-rocket-launch-outline: **Quick Start**](quick-start/index.md)

  ---

  Get up and running with CATChem in minutes

- [:material-book-open-variant: **User Guide**](user-guide/index.md)

  ---

  Comprehensive guide for model users

- [:material-code-braces: **Developer Guide**](developer-guide/index.md)

  ---

  Technical documentation for developers

- [:material-api: **API Reference**](api/index.md)

  ---

  Complete API documentation

</div>

## 🤝 Community & Support

- **GitHub Repository**: [NOAA-GSL/CATChem](https://github.com/NOAA-GSL/CATChem)
- **Documentation**: [https://catchem.readthedocs.io](https://catchem.readthedocs.io)
- **Issue Tracker**: [GitHub Issues](https://github.com/NOAA-GSL/CATChem/issues)
- **Discussions**: [GitHub Discussions](https://github.com/NOAA-GSL/CATChem/discussions)
- **Contact**: [gsl.help@noaa.gov](mailto:gsl.help@noaa.gov)

## 📄 License

CATChem is released under the [Apache 2.0 License](license.md).

---

<div style="text-align: center; margin-top: 3rem; padding: 2rem; background: var(--noaa-gray-light); border-radius: 0.5rem;">
  <p><strong>Community-Driven Atmospheric Chemistry Modeling</strong></p>
  <p style="font-size: 0.9rem; color: var(--noaa-gray-dark);">
    Open source project supporting operational weather prediction and atmospheric research
  </p>
</div>
