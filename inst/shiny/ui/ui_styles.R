# CSS styles and JavaScript for the PerturbPlan dashboard

create_styles <- function() {
  list(
    useShinyjs(),
    
    # Add custom CSS to match shinyapps.io design
    tags$style(HTML("
      /* Custom color scheme - cohesive palette */
      :root {
        --navbar-color: #5B90BF;
        --sidebar-color: #344A5D;
        --sidebar-width: 350px;
        --primary-color: #5B90BF;
        --success-color: #5BA672;
        --info-color: #6BA8CD;
        --warning-color: #D4A945;
        --danger-color: #C85A5A;
        --light-bg: #F8F9FA;
        --section-header: #4A6B82;
        --border-light: #E3E6EA;
        --text-muted: #6C757D;
      }
      
      /* Top navbar styling to match shinyapps.io */
      .main-header .navbar {
        background-color: var(--navbar-color) !important;
        border-bottom: none !important;
        min-height: 50px !important;
      }
      
      .main-header .navbar-custom-menu > .navbar-nav > li > .dropdown-menu {
        background-color: var(--navbar-color) !important;
      }
      
      /* PerturbPlan title styling to match navbar color and sidebar width */
      .main-header .logo {
        background-color: var(--navbar-color) !important;
        color: #fff !important;
        border-bottom: none !important;
        border-right: 1px solid rgba(255,255,255,0.1) !important;
        width: var(--sidebar-width) !important;
        text-align: center !important;
        font-weight: bold !important;
        font-size: 18px !important;
        letter-spacing: 0.5px !important;
      }
      
      .main-header .logo:hover {
        background-color: var(--navbar-color) !important;
      }
      
      /* Sidebar styling to match shinyapps.io */
      .main-sidebar {
        width: var(--sidebar-width) !important;
        background-color: var(--sidebar-color) !important;
        transition: width 0.3s ease-in-out;
        min-height: 100vh !important;
        height: 100vh !important;
      }
      
      .main-sidebar .sidebar {
        background-color: var(--sidebar-color) !important;
        min-height: 100vh !important;
        height: 100vh !important;
        overflow-y: auto !important;
      }
      
      /* Content wrapper adjustments */
      .content-wrapper, .right-side {
        margin-left: var(--sidebar-width) !important;
        transition: margin-left 0.3s ease-in-out;
        min-height: 100vh !important;
      }
      
      /* Download button positioning - push to right */
      .main-header .navbar-custom-menu {
        float: right !important;
        margin-right: 20px !important;
      }
      
      .main-header .navbar .navbar-custom-menu > .navbar-nav {
        margin: 0 !important;
        float: right !important;
      }
      
      .main-header .navbar .navbar-custom-menu > .navbar-nav > li {
        float: right !important;
      }
      
      /* Download button styling - match Plan button exactly (btn-success) */
      #download_results {
        background-color: var(--success-color) !important;
        border-color: var(--success-color) !important;
        color: white !important;
        border-radius: 4px !important;
        padding: 6px 12px !important;
        margin-top: 8px !important;
        font-weight: normal !important;
        font-size: 14px !important;
      }
      
      #download_results:hover, #download_results:focus, #download_results:active {
        background-color: #4A8C5F !important;
        border-color: #4A8C5F !important;
        color: white !important;
      }
      
      /* When sidebar is collapsed */
      .sidebar-collapse .main-sidebar {
        width: 50px !important;
        background-color: var(--sidebar-color) !important;
      }
      
      .sidebar-collapse .main-sidebar .sidebar {
        background-color: var(--sidebar-color) !important;
        overflow: hidden !important;
      }
      
      .sidebar-collapse .content-wrapper, 
      .sidebar-collapse .right-side {
        margin-left: 50px !important;
      }
      
      .sidebar-collapse .main-header .logo {
        width: 50px !important;
      }
      
      /* Ensure sidebar content is hidden when collapsed */
      .sidebar-collapse .main-sidebar .sidebar-menu > li > a > .sidebar-menu-text {
        display: none;
      }
      
      /* Hide our custom parameter panels when collapsed */
      .sidebar-collapse .main-sidebar .sidebar > div {
        display: none !important;
      }
      
      /* Force sidebar color when collapsed */
      .sidebar-collapse .main-sidebar,
      .sidebar-collapse .main-sidebar *,
      .sidebar-collapse .main-sidebar .sidebar,
      .sidebar-collapse .main-sidebar .sidebar-menu,
      .sidebar-collapse .main-sidebar .sidebar > div {
        background: var(--sidebar-color) !important;
        background-color: var(--sidebar-color) !important;
      }
      
      /* Cover everything with sidebar color overlay */
      .sidebar-collapse .main-sidebar::after {
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        width: 50px !important;
        height: 100vh !important;
        background: var(--sidebar-color) !important;
        z-index: 9999 !important;
        pointer-events: none;
      }
      
      /* Responsive design */
      @media (max-width: 767px) {
        .content-wrapper, .right-side {
          margin-left: 0 !important;
        }
        
        .sidebar-collapse .content-wrapper, 
        .sidebar-collapse .right-side {
          margin-left: 0 !important;
        }
        
        .main-header .logo {
          width: 350px !important;
        }
      }
      
      /* Hide the hamburger menu to prevent collapse issues */
      .main-header .sidebar-toggle {
        display: none !important;
      }
      
      /* Ensure proper viewport height for large desktop monitors */
      html, body {
        height: 100% !important;
        min-height: 100vh !important;
      }
      
      .wrapper {
        min-height: 100vh !important;
        height: 100% !important;
      }
      
      /* ========== COHESIVE UI ELEMENT STYLING ========== */
      
      /* Button color scheme */
      .btn-primary {
        background-color: var(--primary-color) !important;
        border-color: var(--primary-color) !important;
        color: white !important;
      }
      
      .btn-primary:hover, .btn-primary:focus, .btn-primary:active {
        background-color: #4A7CA3 !important;
        border-color: #4A7CA3 !important;
        color: white !important;
      }
      
      .btn-success {
        background-color: var(--success-color) !important;
        border-color: var(--success-color) !important;
        color: white !important;
      }
      
      .btn-success:hover, .btn-success:focus, .btn-success:active {
        background-color: #4A8C5F !important;
        border-color: #4A8C5F !important;
        color: white !important;
      }
      
      .btn-info {
        background-color: var(--info-color) !important;
        border-color: var(--info-color) !important;
        color: white !important;
      }
      
      .btn-info:hover, .btn-info:focus, .btn-info:active {
        background-color: #5A94B3 !important;
        border-color: #5A94B3 !important;
        color: white !important;
      }
      
      .btn-warning {
        background-color: var(--warning-color) !important;
        border-color: var(--warning-color) !important;
        color: white !important;
      }
      
      .btn-warning:hover, .btn-warning:focus, .btn-warning:active {
        background-color: #B8923A !important;
        border-color: #B8923A !important;
        color: white !important;
      }
      
      .btn-danger {
        background-color: var(--danger-color) !important;
        border-color: var(--danger-color) !important;
        color: white !important;
      }
      
      .btn-danger:hover, .btn-danger:focus, .btn-danger:active {
        background-color: #A84848 !important;
        border-color: #A84848 !important;
        color: white !important;
      }
      
      /* Clear/secondary buttons - enhanced visibility */
      .btn-default {
        background-color: #ffffff !important;
        border-color: #CED4DA !important;
        color: #495057 !important;
        font-weight: 500 !important;
      }
      
      .btn-default:hover, .btn-default:focus, .btn-default:active {
        background-color: #E9ECEF !important;
        border-color: #ADB5BD !important;
        color: #212529 !important;
      }
      
      /* Go and Clear buttons - match Plan button style */
      #go_overall, #go_slice {
        background-color: var(--primary-color) !important;
        border-color: var(--primary-color) !important;
        color: white !important;
        font-weight: normal !important;
        font-size: 14px !important;
      }
      
      #go_overall:hover, #go_overall:focus, #go_overall:active,
      #go_slice:hover, #go_slice:focus, #go_slice:active {
        background-color: #4A7CA3 !important;
        border-color: #4A7CA3 !important;
        color: white !important;
      }
      
      /* Clear buttons - use default secondary style but enhanced */
      #clear, #slice_clear {
        background-color: #ffffff !important;
        border-color: #CED4DA !important;
        color: #495057 !important;
        font-weight: normal !important;
        font-size: 14px !important;
      }
      
      #clear:hover, #clear:focus, #clear:active,
      #slice_clear:hover, #slice_clear:focus, #slice_clear:active {
        background-color: #E9ECEF !important;
        border-color: #ADB5BD !important;
        color: #212529 !important;
      }
      
      /* Box/panel styling - standardized to primary blue */
      .box.box-primary .box-header,
      .box.box-info .box-header,
      .box.box-success .box-header,
      .box.box-warning .box-header,
      .box.box-danger .box-header {
        background-color: var(--primary-color) !important;
        color: white !important;
      }
      
      .box.box-primary,
      .box.box-info,
      .box.box-success,
      .box.box-warning,
      .box.box-danger {
        border-top-color: var(--primary-color) !important;
      }
      
      /* Generic box header styling for consistency */
      .box-header {
        background-color: var(--primary-color) !important;
        color: white !important;
      }
      
      .box {
        border-top-color: var(--primary-color) !important;
      }
      
      /* Collapsible section headers in sidebar */
      #exp-header, #perturbation-header, #analysis-header, #effects-header {
        background-color: var(--section-header) !important;
        color: white !important;
        border-bottom: 1px solid rgba(255,255,255,0.2) !important;
      }
      
      #exp-header:hover, #perturbation-header:hover, #analysis-header:hover, #effects-header:hover {
        background-color: #3E5A73 !important;
      }
      
      #exp-chevron, #perturbation-chevron, #analysis-chevron, #effects-chevron {
        color: white !important;
      }
      
      /* Section containers - reduced spacing */
      #experimental-content, #perturbation-content, #analysis-content, #effects-content {
        background-color: var(--light-bg) !important;
        border: 1px solid var(--border-light) !important;
        border-radius: 0 0 4px 4px !important;
        color: #212529 !important;
        padding: 10px !important;
      }
      
      /* Ensure all text within sidebar sections is dark and readable */
      #experimental-content *, #perturbation-content *, #analysis-content *, #effects-content * {
        color: #212529 !important;
      }
      
      #experimental-content label, #perturbation-content label, #analysis-content label, #effects-content label {
        color: #212529 !important;
        font-weight: 500 !important;
      }
      
      /* Reduce vertical spacing in sidebar form elements */
      .main-sidebar .form-group {
        margin-bottom: 10px !important;
      }
      
      /* Reduce margin between form inputs */
      .main-sidebar .form-control {
        margin-bottom: 5px !important;
      }
      
      /* Reduce spacing between collapsible sections */
      .main-sidebar > div > div[style*='margin-bottom: 10px'] {
        margin-bottom: 5px !important;
      }
      
      /* Reduce spacing in conditionalPanel elements */
      .main-sidebar .shiny-input-container {
        margin-bottom: 8px !important;
      }
      
      /* Horizontal separator line in sidebar */
      .sidebar-separator {
        border: none !important;
        border-top: 1px solid rgba(255, 255, 255, 0.2) !important;
        margin: 15px 0 !important;
        background: none !important;
      }
      
      /* Center Plan button horizontally in sidebar */
      #plan_btn {
        display: block !important;
        margin: 0 auto !important;
      }
      
      /* Form inputs and labels styling */
      .form-control:focus {
        border-color: var(--primary-color) !important;
        box-shadow: 0 0 0 0.2rem rgba(91, 144, 191, 0.25) !important;
      }
      
      .form-control {
        border-color: var(--border-light) !important;
      }
      
      /* Ensure form labels are visible with dark text */
      .control-label, label, .checkbox label, .radio label {
        color: #212529 !important;
        font-weight: 500 !important;
      }
      
      /* Sidebar form labels specifically */
      .main-sidebar label, .main-sidebar .control-label {
        color: #212529 !important;
        font-weight: 500 !important;
      }
      
      /* Form group text */
      .form-group label {
        color: #212529 !important;
        font-weight: 500 !important;
      }
      
      /* Help text and descriptions */
      .help-block, .form-text {
        color: var(--text-muted) !important;
      }
      
      /* Radio buttons and checkboxes */
      .radio input[type='radio']:checked + label::before,
      .checkbox input[type='checkbox']:checked + label::before {
        background-color: var(--primary-color) !important;
        border-color: var(--primary-color) !important;
      }
      
      /* Alert/banner styling */
      .alert-info {
        background-color: rgba(107, 168, 205, 0.1) !important;
        border-color: var(--info-color) !important;
        color: #2E5C73 !important;
      }
      
      .alert-success {
        background-color: rgba(91, 166, 114, 0.1) !important;
        border-color: var(--success-color) !important;
        color: #2E5C3A !important;
      }
      
      .alert-warning {
        background-color: rgba(212, 169, 69, 0.1) !important;
        border-color: var(--warning-color) !important;
        color: #6B5425 !important;
      }
      
      .alert-danger {
        background-color: rgba(200, 90, 90, 0.1) !important;
        border-color: var(--danger-color) !important;
        color: #6B2E2E !important;
      }
      
      /* Tab styling */
      .nav-tabs > li.active > a, .nav-tabs > li.active > a:focus, .nav-tabs > li.active > a:hover {
        background-color: var(--light-bg) !important;
        border-bottom-color: var(--primary-color) !important;
        color: var(--primary-color) !important;
      }
      
      .nav-tabs > li > a:hover {
        border-color: var(--border-light) var(--border-light) var(--primary-color) !important;
        color: var(--primary-color) !important;
      }
      
      /* Custom file upload info banner */
      .file-upload-info {
        background-color: rgba(212, 169, 69, 0.1) !important;
        border: 1px solid var(--warning-color) !important;
        color: #6B5425 !important;
      }
      
      .file-upload-success {
        background-color: rgba(91, 166, 114, 0.1) !important;
        border: 1px solid var(--success-color) !important;
        color: #2E5C3A !important;
      }
    ")),
    
    # Add custom JavaScript for collapsible sections
    tags$script(HTML("
      function toggleSection(contentId, chevronId) {
        console.log('toggleSection called with:', contentId, chevronId);
        
        // First, collapse all other sections
        var allSections = ['experimental-content', 'perturbation-content', 'analysis-content', 'effects-content'];
        var allChevrons = ['exp-chevron', 'perturbation-chevron', 'analysis-chevron', 'effects-chevron'];
        
        for (var i = 0; i < allSections.length; i++) {
          var section = document.getElementById(allSections[i]);
          var chevron = document.getElementById(allChevrons[i]);
          
          if (allSections[i] === contentId) {
            // Toggle the clicked section
            if (section.style.display === 'none') {
              section.style.display = 'block';
              chevron.className = 'fa fa-chevron-down';
            } else {
              section.style.display = 'none';
              chevron.className = 'fa fa-chevron-right';
            }
          } else {
            // Collapse other sections
            section.style.display = 'none';
            chevron.className = 'fa fa-chevron-right';
          }
        }
      }
      
      // Initialize on page load
      $(document).ready(function() {
        console.log('JavaScript loaded');
        
        // Set initial states explicitly
        setTimeout(function() {
          var exp = document.getElementById('experimental-content');
          var perturbation = document.getElementById('perturbation-content');
          var analysis = document.getElementById('analysis-content');
          var effects = document.getElementById('effects-content');
          
          var expChevron = document.getElementById('exp-chevron');
          var perturbationChevron = document.getElementById('perturbation-chevron');
          var analysisChevron = document.getElementById('analysis-chevron');
          var effectsChevron = document.getElementById('effects-chevron');
          
          console.log('Elements found:', exp ? 'exp-yes' : 'exp-no', 
                                       perturbation ? 'perturbation-yes' : 'perturbation-no',
                                       analysis ? 'analysis-yes' : 'analysis-no',
                                       effects ? 'effects-yes' : 'effects-no');
          
          if (exp && perturbation && analysis && effects) {
            // Set initial state: only experimental section open
            exp.style.display = 'block';
            perturbation.style.display = 'none';
            analysis.style.display = 'none';
            effects.style.display = 'none';
            
            // Set chevron states
            if (expChevron) expChevron.className = 'fa fa-chevron-down';
            if (perturbationChevron) perturbationChevron.className = 'fa fa-chevron-right';
            if (analysisChevron) analysisChevron.className = 'fa fa-chevron-right';
            if (effectsChevron) effectsChevron.className = 'fa fa-chevron-right';
            
            console.log('Initial states set');
          } else {
            console.log('Could not find all elements');
          }
        }, 1000); // Wait 1 second for Shiny to fully load
      });
    "))
  )
}