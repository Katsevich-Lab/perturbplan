# CSS styles and JavaScript for the PerturbPlan dashboard

create_styles <- function() {
  list(
    useShinyjs(),
    
    # Add custom CSS to adjust sidebar width and color with proper collapse support
    tags$style(HTML("
      .main-sidebar {
        width: 350px !important;
        background-color: #3c8dbc !important;
        transition: width 0.3s ease-in-out;
      }
      
      .main-sidebar .sidebar {
        background-color: #3c8dbc !important;
      }
      
      .content-wrapper, .right-side {
        margin-left: 350px !important;
        transition: margin-left 0.3s ease-in-out;
      }
      
      /* When sidebar is collapsed */
      .sidebar-collapse .main-sidebar {
        width: 50px !important;
        background-color: #3c8dbc !important;
      }
      
      .sidebar-collapse .main-sidebar .sidebar {
        background-color: #3c8dbc !important;
        overflow: hidden !important;
      }
      
      .sidebar-collapse .content-wrapper, 
      .sidebar-collapse .right-side {
        margin-left: 50px !important;
      }
      
      /* Ensure sidebar content is hidden when collapsed */
      .sidebar-collapse .main-sidebar .sidebar-menu > li > a > .sidebar-menu-text {
        display: none;
      }
      
      /* Hide our custom parameter panels when collapsed */
      .sidebar-collapse .main-sidebar .sidebar > div {
        display: none !important;
      }
      
      /* Nuclear option - force blue background everywhere */
      .sidebar-collapse .main-sidebar,
      .sidebar-collapse .main-sidebar *,
      .sidebar-collapse .main-sidebar .sidebar,
      .sidebar-collapse .main-sidebar .sidebar-menu,
      .sidebar-collapse .main-sidebar .sidebar > div {
        background: #3c8dbc !important;
        background-color: #3c8dbc !important;
      }
      
      /* Cover everything with blue overlay */
      .sidebar-collapse .main-sidebar::after {
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        width: 50px !important;
        height: 100vh !important;
        background: #3c8dbc !important;
        z-index: 9999 !important;
        pointer-events: none;
      }
      
      @media (max-width: 767px) {
        .content-wrapper, .right-side {
          margin-left: 0 !important;
        }
        
        .sidebar-collapse .content-wrapper, 
        .sidebar-collapse .right-side {
          margin-left: 0 !important;
        }
      }
      
      /* Hide the hamburger menu to prevent collapse issues */
      .main-header .sidebar-toggle {
        display: none !important;
      }
    ")),
    
    # Add custom JavaScript for collapsible sections
    tags$script(HTML("
      function toggleSection(contentId, chevronId) {
        console.log('toggleSection called with:', contentId, chevronId);
        
        // First, collapse all other sections
        var allSections = ['experimental-content', 'analysis-content', 'effects-content'];
        var allChevrons = ['exp-chevron', 'analysis-chevron', 'effects-chevron'];
        
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
          var analysis = document.getElementById('analysis-content');
          var effects = document.getElementById('effects-content');
          
          var expChevron = document.getElementById('exp-chevron');
          var analysisChevron = document.getElementById('analysis-chevron');
          var effectsChevron = document.getElementById('effects-chevron');
          
          console.log('Elements found:', exp ? 'exp-yes' : 'exp-no', 
                                       analysis ? 'analysis-yes' : 'analysis-no',
                                       effects ? 'effects-yes' : 'effects-no');
          
          if (exp && analysis && effects) {
            // Set initial state: only experimental section open
            exp.style.display = 'block';
            analysis.style.display = 'none';
            effects.style.display = 'none';
            
            // Set chevron states
            if (expChevron) expChevron.className = 'fa fa-chevron-down';
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