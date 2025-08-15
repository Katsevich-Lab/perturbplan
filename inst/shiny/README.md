# PerturbPlan Shiny App - Modular Structure

This directory contains the modularized Shiny application for PerturbPlan power analysis.

## File Structure

```
inst/shiny/
├── app.R                 # Main entry point
├── ui/
│   ├── ui_main.R        # Main UI composition
│   ├── ui_header.R      # Dashboard header with download button
│   ├── ui_sidebar.R     # Parameter input sidebar with collapsible sections
│   ├── ui_styles.R      # CSS styles and JavaScript
│   └── ui_tabs.R        # Tab panels for different views
└── server/
    ├── server_main.R    # Main server composition  
    ├── server_power.R   # Power calculation logic
    ├── server_selection.R # Selection and interaction handling
    ├── server_plots.R   # Plot rendering (heatmap, slice)
    └── server_curves.R  # Power curves and download functionality
```

## UI Modules

### `ui_header.R`
- Dashboard header with "PerturbPlan" title
- Conditional download button that appears after planning

### `ui_sidebar.R` 
- Collapsible parameter sections:
  - Experimental choices (MOI, targets, gRNAs, etc.)
  - Analysis choices (TPM threshold, FDR, test side, etc.)
  - Assumed effect sizes (fold-change parameters)
- Centered "Plan" button

### `ui_styles.R`
- Custom CSS for sidebar width (350px) and blue color theme
- Hamburger menu disabled to prevent collapse issues  
- JavaScript for collapsible section functionality
- Auto-collapse logic for completed parameter sections

### `ui_tabs.R`
- Three main tabs in content area:
  - "Overall power (heatmap)" - Interactive heatmap view
  - "Overall power (slice)" - Line plots for drilling down
  - "Power over TPM & FC" - Detailed power curves

## Server Modules

### `server_power.R`
- Core power calculation using `perturbplan::calculate_power_grid()`
- Reactive values for planned state, power results, grid data
- Helper functions for selection logic

### `server_selection.R`
- Handles user interactions with plots (clicks, selections)
- Manages selection state for rows, columns, and tiles
- UI generation for selection controls

### `server_plots.R`
- Renders heatmap with selection rectangles
- Slice plot generation and interactions
- Multi-selection support for drill-down analysis

### `server_curves.R`
- On-demand power curve computation for selected tiles
- Combined TPM and fold-change curve plots using patchwork
- Excel download functionality with multiple sheets

## Key Features

1. **Modular Design**: Each component is self-contained and focused
2. **Data Flow**: Modules communicate through returned reactive values
3. **Reusability**: Components can be easily modified or extended
4. **Maintainability**: Clear separation of concerns
5. **Testing**: Each module can be tested independently

## Usage

The modular structure maintains all existing functionality while improving code organization:

- Launch with `perturbplan::launch_app()`
- All parameter inputs, visualizations, and download features work as before
- Code is now easier to maintain and extend