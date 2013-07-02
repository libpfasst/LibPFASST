#
# LIBPFASST documentation build configuration file
#

import sys, os

# path to autogen'ed modules
sys.path.append(os.path.abspath('.'))
sys.path.append(os.path.abspath('..'))

# extentions
extensions = ['sphinx.ext.autodoc', 'mathjax', 'sphinx.ext.coverage']
mathjax_path = 'http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_HTML'

# general configuration
source_suffix  = '.rst'
master_doc     = 'index'

# html configuration
pygments_style = 'sphinx'
html_theme     = 'default'
html_theme_options = {
    'stickysidebar': True,
    }
#html_short_title = 'PFASST'
html_show_sourcelink = False
html_domain_indices = False

html_sidebars = {
    '**': ['globaltoc.html', 'searchbox.html'],
    }

# project information
project   = 'LIBPFASST'
copyright = '2011, 2012, 2013 Matthew Emmett'

release = 0.5.0
