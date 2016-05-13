
# Written by Chris Calderon 2015 (Chris.Calderon@UrsaAnalytics.com)
#
#
# Copyright 2015 Ursa Analytics, Inc.

   # Licensed under the Apache License, Version 2.0 (the "License");
   # you may not use this file except in compliance with the License.
   # You may obtain a copy of the License at

   #     http://www.apache.org/licenses/LICENSE-2.0
   
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


__license__ = "Apache License, Version 2.0"
__author__  = "Chris Calderon, Ursa Analytics, Inc. [www.UrsaAnalytics.com]"
__status__  = "Development"


"""
Ursa_IPyNBpltWrapper.py

Module containing some functions allowing easier matplotlib plotting in IPyNB [hides some messy details in matplotlib func calls]
Also contains some I/O routines commonly used when analyzing data.

Written by Chris Calderon 2015 (Chris.Calderon@UrsaAnalytics.com)


Copyright 2015 Ursa Analytics, Inc.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

"""

import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np


defaultFont = 'Bitstream Vera Sans'
# defaultFont = 'Helvetica' #this font often requires mods to system to allow python to find font.  distribute code with a "prepackaged" font

# def basicPltSetup(axisfontsize='22',legendfontsize='18',pltFont='Helvetica'):  #Helvetica font isn't bundled with most distros so requires installation if desired 
def basicPltSetup(axisfontsize='22',legendfontsize='18',pltFont=defaultFont):
    
    """
    Simple wrapper for using custom font (e.g., Helvetica)
    and desired axis tick / font sizes (returns a dictionary)

    For distribution, probably want to remove fontname since users might have a hard time installing custom fonts 
    """
    pltPars={}
    axis_font = {'fontname':pltFont, 'size':axisfontsize}
    title_font = {'fontname':pltFont, 'size':axisfontsize, 'color':'black', 'weight':'normal',
                  'verticalalignment':'bottom'} # Bottom vertical alignment for more space

    #commands below don't easily port to other operating systems.  just comment out features using this tweak.
    # font_path = '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/' + pltFont +'.ttf'
    # font_prop = font_manager.FontProperties(fname=font_path, size=axisfontsize) #used in custom legend below 
    # legendfont_prop = font_manager.FontProperties(fname=font_path, size=legendfontsize) #used in custom legend below
    
    fig = plt.figure() 
    ax  = fig.add_subplot(111,adjustable='box')

    #pack things up for revising plot later in dict (see if possible to adjust before plot made)
    pltPars['axis_font'] = axis_font

    # pltPars['leg_prop'] = legendfont_prop #e.g. of use: lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,ncol=1,numpoints=1,prop=pltPars['leg_prop'])
    pltPars['axisfontsize'] = axisfontsize

    return fig,ax,pltPars

# def sharedXPltSetup(axisfontsize='22',legendfontsize='18',pltFont='Helvetica'):  #Helvetica font isn't bundled with most distros so requires installation if desired 
def sharedXPltSetup(axisfontsize='22',legendfontsize='18',pltFont=defaultFont):  

    """
    Simple wrapper for setting up shared x-axis with custom font (e.g., Helvetica) and desired axis tick / font sizes (returns a dictionary)
    After adding data to plots, likely need to rerun
    fig.subplots_adjust(hspace=0) 
    in ipythonNB


    For distribution, probably want to remove fontname since users might have a hard time installing custom fonts 
    """
    
    pltPars={}
    axis_font = {'fontname':pltFont, 'size':axisfontsize}
    title_font = {'fontname':pltFont, 'size':axisfontsize, 'color':'black', 'weight':'normal',
                  'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    
    # font_path = '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/' + pltFont +'.ttf'
    # font_prop = font_manager.FontProperties(fname=font_path, size=axisfontsize) #used in custom legend below 
    # legendfont_prop = font_manager.FontProperties(fname=font_path, size=legendfontsize) #used in custom legend below

    fig = plt.figure() 
    fig, (axTop, ax) = plt.subplots(2, sharex=True)
    fig.subplots_adjust(hspace=0) #cannot think of situation where this isn't desirable in shared x-axis setting

    #pack things up for revising plot later in dict (see if possible to adjust before plot made)
    pltPars['axis_font'] = axis_font
    pltPars['title_font'] = title_font
    # pltPars['font_prop'] = font_prop
    # pltPars['leg_prop'] = legendfont_prop #e.g. of use: lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,ncol=1,numpoints=1,prop=pltPars['leg_prop'])
    pltPars['axisfontsize'] = axisfontsize

    return fig,ax,axTop,pltPars

def xyLabels(ax,pltPars,xlabel='',ylabel='',title='',pltFont=defaultFont): #don't use Helvetic in deployed versions (likely a pain to setup)
    """
    set the axis number font size and type.  option, put  x and y labels [for pubs, leave empty and use latex overlay ] 
    (requires pltPar dict from basicPltSetup) 
    """
    ax.set_xlabel(xlabel, **pltPars['axis_font'])
    ax.set_ylabel(ylabel, **pltPars['axis_font'])
    ax.set_title(title, **pltPars['axis_font'])
    
    #set font size
    for item in ([ax.xaxis.label, ax.yaxis.label ] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(pltPars['axisfontsize'])
    
    # Set the tick font type
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontname(pltFont)




def wout(floc,DATAMAT): #write out an ASCII file to "floc" given a 2D numpy array "DATAMAT"
    """
    write an ASCII file.  

    floc: file location  (can use "~" in path) 
    DATAMAT: numpy matrix or vector
    """
    if floc[0]=='~':  #replace standard unix home shortcut with explicit path
                floc='/Users/calderoc' + floc[1:]
    np.savetxt(floc, DATAMAT)
    return 0

def randomColors(nsamps):
    """
    return of list of length nsamps containing randomly selected colors for line plots.  nice alternative to boring
    #defaults (useful for plotting many lines on one graph)
    # b: blue
    # g: green
    # r: red
    # c: cyan
    # m: magenta
    # y: yellow
    # k: black
    # w: white

    matplotlib.colors.cnames allows things like "orange" and "fuchsia"

    """ 
    clist = matplotlib.colors.cnames #create master list of default plotting colors
    myclrs = [] #create a list storing a random selection of default color names stored in matplotlib.colors.cnames
    for i in range(nsamps):
        myclrs.append(random.choice(list(clist.keys()))) #pick a random key value
    return myclrs



