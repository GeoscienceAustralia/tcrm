#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Title: configeditor.py
Author: Nicholas Summons, nicholas.summons@ga.gov.au
CreationDate: 2010-09-09
Description: A graphical user interface for generating and editing TCRM configuration files.
"""
import os
import wx
import time
import wx.richtext
import wx.grid
import csv
import math
import numpy
import unicodedata
import ConfigParser
from copy import deepcopy
from wx import ImageFromStream, BitmapFromImage, EmptyIcon
import cStringIO, zlib
from Utilities import pathLocator
import sqlite3

# Switch off minor warning messages
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="pytz")


class MainFrame ( wx.Frame ):

    def __init__( self, parent ):

        tcrm_dir = pathLocator.getRootDirectory()
        tcrm_input_dir = os.path.join(tcrm_dir, 'input')
        localitiesDataFile = os.path.join(tcrm_input_dir, 'localities.dat')
 
        self.sqlcon = sqlite3.connect(localitiesDataFile)
        self.sqlcur = self.sqlcon.cursor()

        # Sort country names into alphabetical order (since dictionaries do not preserve order)
        self.sqlcur.execute('select placename from localities where placetype=?', ('country',))
        allStations_keys = [z[0] for z in self.sqlcur.fetchall()]
        allStations_keys.sort()

        self.columnChoices = {'refName':['', 'num', 'date', 'season', 'index', 'year', 'month', 'day',
                                         'hour', 'minute', 'lat', 'lon', 'pressure', 'rmax', 'tcserialno'],
                              'displayStr':['', 'Storm Number', 'Date', 'Season', 'Index', 'Year', 'Month', 'Day',
                                            'Hour', 'Minute', 'Latitude', 'Longitude', 'Central Pressure', 'RMW', 'TC Serial Number'],
                              'unitMapping':['', '', 'DateFormat', '', '', '', '', '',
                                             '', '', '', '', 'PressureUnits', 'LengthUnits', ''],
                              'Countries':[u''] + allStations_keys,
                              'Divisions':'',
                              'Locations':''}
        self.speedUnits = ['mps', 'mph', 'kph', 'kts']
        self.speedUnitsLongname = ['metres / second','miles / hour','kilometres / hour','knots']
        self.pressureUnits = ['hPa', 'kPa', 'Pa', 'inHg', 'mmHg']
        self.lengthUnits = ['km', 'mi', 'nm']
        self.KDETypes = ['Biweight', 'Epanechnikov', 'Triangular', 'Gaussian']
        self.profileTypes = ['powell', 'holland']
        self.profileTypes_display = ['Powell', 'Holland']
        self.windFieldTypes = ['kepert']
        self.windFieldTypes_display = ['Kepert']
        self.tcrmPath = pathLocator.getRootDirectory()
        self.tcrmInputDir = os.path.join(self.tcrmPath, 'input')
        initialmsg = " Welcome to the Tropical Cyclone Risk Model (TCRM) configuration editor.\n" + \
                     " To start, press the 'Create New' button to create a new configuration file or\n" + \
                     " 'Browse' to open a pre-existing file."
        self.userhelptext = []
        self.gridMaxRows = 20
        self.gridRowSize = 17
        self.gridMaxCols = 18
        self.gridColSize = 100
        self.datafilename = ''
        self.defaultconfigname = 'untitled.ini'
        self.defaultSourceName = 'UnnamedDataSource'
        self.configDict = {}
        self.configDictSections = []
        self.nextRecNo = 0
        self.tabsEnabled = False
        self.settingsChanged = False
        self.valueWarningRaised = False

        self.display_stats_tab = False
        self.display_RMW_tab = False
        self.display_MSLP_tab = False
        self.display_windfield_tab = False

        wx.Frame.__init__  ( self, parent, id = wx.ID_ANY, title = u"Tropical Cyclone Risk Model - Configuration Editor", pos = wx.DefaultPosition, size = wx.Size( 574,630 ),
                             style=wx.MINIMIZE_BOX|wx.SYSTEM_MENU|wx.CAPTION|wx.TAB_TRAVERSAL|wx.CLOSE_BOX|wx.CLIP_CHILDREN)

        self.Bind(wx.EVT_CLOSE, self.OnProjectExit)

        self.pointers = []

        icon = getIcon()
        self.SetIcon(icon)

        self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )

        self.MenuBar = wx.MenuBar( 0 )
        self.File = wx.Menu()
        self.file_new = wx.MenuItem( self.File, wx.ID_ANY, u"&New", u" Create a new configuration file", wx.ITEM_NORMAL )
        self.File.AppendItem( self.file_new )
        self.Bind(wx.EVT_MENU, self.menu_sel_file_new, self.file_new)

        self.file_open = wx.MenuItem( self.File, wx.ID_ANY, u"Open", u" Open an existing configuration file", wx.ITEM_NORMAL )
        self.File.AppendItem( self.file_open )
        self.Bind(wx.EVT_MENU, self.menu_sel_file_open, self.file_open)

        self.file_save = wx.MenuItem( self.File, wx.ID_ANY, u"Save", u" Save configuration file", wx.ITEM_NORMAL)
        self.File.AppendItem( self.file_save )
        self.file_save.Enable( False )
        self.Bind(wx.EVT_MENU, self.menu_sel_file_save, self.file_save)

        self.file_saveas = wx.MenuItem( self.File, wx.ID_ANY, u"Save As", u" Save configuration file as ......", wx.ITEM_NORMAL )
        self.File.AppendItem( self.file_saveas )
        self.file_saveas.Enable( False )
        self.Bind(wx.EVT_MENU, self.menu_sel_file_saveas, self.file_saveas)

        self.File.AppendSeparator()

        self.file_exit = wx.MenuItem( self.File, wx.ID_ANY, u"Exit", u" Exit TCRM", wx.ITEM_NORMAL )
        self.File.AppendItem( self.file_exit )
        self.Bind(wx.EVT_MENU, self.menu_sel_file_exit, self.file_exit)

        self.MenuBar.Append( self.File, u"File" )

        self.Help = wx.Menu()

        self.help_about = wx.MenuItem( self.Help, 1010, u"About", " About TCRM", wx.ITEM_NORMAL )
        self.Help.AppendItem( self.help_about )
        self.Bind(wx.EVT_MENU, self.loadHelp, self.help_about)

        self.MenuBar.Append( self.Help, u"Help" )

        self.SetMenuBar( self.MenuBar )

        sizer00 = wx.BoxSizer( wx.VERTICAL )

        self.backpanel = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.guiFont = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)
        self.guiFont.SetPointSize(8)
        self.backpanel.SetFont(self.guiFont)
        sizer0 = wx.BoxSizer( wx.VERTICAL )

        bSizer42 = wx.BoxSizer( wx.VERTICAL )
        sizer0.Add( bSizer42, 1, wx.EXPAND, 5 )

        sizer01 = wx.BoxSizer( wx.HORIZONTAL )

        self.label011 = wx.StaticText( self.backpanel, wx.ID_ANY, u"Configuration File:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label011.Wrap( -1 )
        sizer01.Add( self.label011, 0, wx.ALIGN_CENTER_VERTICAL, 5 )

        self.textCtrl012 = wx.TextCtrl( self.backpanel, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_PROCESS_ENTER|wx.TE_READONLY )
        sizer01.Add( self.textCtrl012, 1, wx.ALIGN_CENTER|wx.ALL, 5 )

        self.button013 = wx.Button( self.backpanel, wx.ID_ANY, u"Create New", wx.DefaultPosition, wx.DefaultSize, 0 )
        sizer01.Add( self.button013, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        self.button013.Bind(wx.EVT_BUTTON, self.createNewConfigFile)

        self.button014 = wx.Button( self.backpanel, wx.ID_ANY, u"Browse", wx.DefaultPosition, wx.DefaultSize, 0 )
        sizer01.Add( self.button014, 0, wx.ALIGN_CENTER|wx.ALL, 5 )
        self.button014.Bind(wx.EVT_BUTTON, self.LoadConfigFile)

        sizer0.Add( sizer01, 6, wx.ALL|wx.EXPAND, 13 )

        staticBox02 = wx.StaticBoxSizer( wx.StaticBox( self.backpanel, wx.ID_ANY, u"Settings" ), wx.VERTICAL )

        self.notebook_panel = wx.Notebook( self.backpanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0 )


        #-------------------- Create Input Tab -------------------
        self.tab_input = wx.Panel( self.notebook_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        sizer1 = wx.BoxSizer( wx.VERTICAL )

        staticBox11 = wx.StaticBoxSizer( wx.StaticBox( self.tab_input, wx.ID_ANY, u"Input File" ), wx.HORIZONTAL )

        self.label111 = wx.StaticText( self.tab_input, wx.ID_ANY, u"Input CSV Track File:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label111.Wrap( -1 )
        staticBox11.Add( self.label111, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.textCtrl1119 = wx.TextCtrl( self.tab_input, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_PROCESS_ENTER|wx.TE_READONLY )
        staticBox11.Add( self.textCtrl1119, 1, wx.ALIGN_CENTER|wx.ALL, 5 )

        self.m_button11118 = wx.Button(self.tab_input, wx.ID_ANY, u"Browse", wx.DefaultPosition, wx.DefaultSize, 0)
        staticBox11.Add( self.m_button11118, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        self.m_button11118.Bind(wx.EVT_BUTTON, self.SelectCSVFile)

        sizer1.Add( staticBox11, 1, wx.ALL|wx.EXPAND, 5 )

        staticBox12 = wx.StaticBoxSizer( wx.StaticBox( self.tab_input, wx.ID_ANY, u"Map Columns" ), wx.VERTICAL )

        self.label121 = wx.StaticText( self.tab_input, wx.ID_ANY, u"Select Columns from Drop Down Menus:", wx.Point( -1,-1 ), wx.DefaultSize, 0 )
        self.label121.Wrap( -1 )
        staticBox12.Add( self.label121, 0, wx.ALL, 7 )

        self.scrolledWindow122 = wx.ScrolledWindow( self.tab_input, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.HSCROLL|wx.STATIC_BORDER|wx.VSCROLL )
        self.scrolledWindow122.SetScrollRate( 5, 5 )
        Sizer122 = wx.BoxSizer( wx.VERTICAL )

        Sizer1221 = wx.BoxSizer( wx.HORIZONTAL )

        # Create drop down menus for input file column mapping
        self.colmap = []

        for k in range(self.gridMaxCols):
            choicemenu = wx.Choice(self.scrolledWindow122, 9000 + k, wx.DefaultPosition, wx.DefaultSize, self.columnChoices['displayStr'], 0)
            choicemenu.SetSelection(0)
            choicemenu.Enable(True)
            choicemenu.SetMinSize(wx.Size(self.gridColSize,-1))
            self.Bind(wx.EVT_CHOICE, self.onColMapChange, choicemenu)
            self.colmap.append(choicemenu)
            Sizer1221.Add(self.colmap[k], 0, 0, 5)

        Sizer122.Add( Sizer1221, 0, 0, 5 )

        self.grid_inputCSV = wx.grid.Grid( self.scrolledWindow122, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0 )

        # Grid
        self.grid_inputCSV.CreateGrid( self.gridMaxRows, self.gridMaxCols )
        self.grid_inputCSV.EnableEditing( False )
        self.grid_inputCSV.EnableGridLines( True )
        self.grid_inputCSV.EnableDragGridSize( False )
        self.grid_inputCSV.SetMargins( 0, 0 )

        # Columns
        for k in range(self.gridMaxCols):
            self.grid_inputCSV.SetColSize( k, self.gridColSize )
        self.grid_inputCSV.EnableDragColMove( False )
        self.grid_inputCSV.EnableDragColSize( False )
        self.grid_inputCSV.SetColLabelSize( 0 )
        self.grid_inputCSV.SetColLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )

        # Rows
        for k in range(self.gridMaxRows):
            self.grid_inputCSV.SetRowSize( k, self.gridRowSize )
        self.grid_inputCSV.EnableDragRowSize( False )
        self.grid_inputCSV.SetRowLabelSize( 0 )
        self.grid_inputCSV.SetRowLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )

        # Cell Defaults
        self.grid_inputCSV.SetDefaultCellAlignment( wx.ALIGN_RIGHT, wx.ALIGN_TOP )
        Sizer122.Add( self.grid_inputCSV, 0, 0, 5 )

        self.scrolledWindow122.SetSizer( Sizer122 )
        self.scrolledWindow122.Layout()
        Sizer122.Fit( self.scrolledWindow122 )
        staticBox12.Add( self.scrolledWindow122, 1, wx.ALL|wx.EXPAND, 5 )

        sizerhz1 = wx.BoxSizer( wx.HORIZONTAL )

        sizerhz2 = wx.BoxSizer( wx.HORIZONTAL )
        sizerhz1.Add( sizerhz2, 6, wx.ALL, 7 )

        self.labeln0 = wx.StaticText( self.tab_input, wx.ID_ANY, u"Number of text header lines:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.labeln0.Wrap( -1 )
        sizerhz1.Add( self.labeln0, 0, wx.ALL, 7 )

        self.m_spinCtrl1 = wx.SpinCtrl( self.tab_input, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0, 0, 10, 0 )
        self.m_spinCtrl1.Bind(wx.EVT_SPINCTRL, self.headerLinesSettingChange)
        sizerhz1.Add( self.m_spinCtrl1, 1, wx.ALL, 1 )

        staticBox12.Add( sizerhz1, 0, wx.EXPAND, 5 )

        sizer1.Add( staticBox12, 2, wx.ALL|wx.EXPAND, 5 )

        self.tab_input.SetSizer( sizer1 )
        self.tab_input.Layout()
        sizer1.Fit( self.tab_input )
        self.notebook_panel.AddPage( self.tab_input, u"Input", True )
        self.userhelptext.append(" Please select an input track file by pressing the 'Browse' button located to \n" + \
                                 " the right of 'Input CSV Track File'. " + \
                                 "\n\n                                          { Scroll down for further help }      " + \
                                 "\n\n Once selected, sample data will populate the columns above." + \
                                 "\n Use the drop-down column menus to instruct TCRM what each column" + \
                                 "\n represents.  Any columns not required by TCRM can be left with " + \
                                 "\n blank assignments.")
        #---------------------------------------------------------


        #------------------- Create Output Tab -------------------
        self.tab_output = wx.Panel( self.notebook_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        sizer2 = wx.BoxSizer( wx.VERTICAL )

        staticBox21 = wx.StaticBoxSizer( wx.StaticBox( self.tab_output, wx.ID_ANY, u"Output Location" ), wx.HORIZONTAL )

        self.label211 = wx.StaticText( self.tab_output, wx.ID_ANY, u"Output Directory: ", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label211.Wrap( -1 )
        staticBox21.Add( self.label211, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.textCtrl212 = wx.TextCtrl( self.tab_output, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_PROCESS_ENTER|wx.TE_READONLY )
        staticBox21.Add( self.textCtrl212, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.button213 = wx.Button( self.tab_output, wx.ID_ANY, u"Browse", wx.DefaultPosition, wx.DefaultSize, 0 )
        staticBox21.Add( self.button213, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        self.button213.Bind(wx.EVT_BUTTON, self.LoadDirPicker)
        
        sizer2.Add( staticBox21, 1, wx.ALL|wx.EXPAND, 5 )

        padder22 = wx.BoxSizer( wx.VERTICAL )

        sizer2.Add( padder22, 2, wx.EXPAND, 5 )

        self.tab_output.SetSizer( sizer2 )
        self.tab_output.Layout()
        sizer2.Fit( self.tab_output )
        self.notebook_panel.AddPage( self.tab_output, u"Output", False )
        self.userhelptext.append(" Press the 'Browse' button to select a directory for the output data.")
        #---------------------------------------------------------


        #------------------- Create Domain Tab -------------------
        self.tab_domain = wx.Panel( self.notebook_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.tab_domain.SetMinSize( wx.Size( -1,1 ) )

        sizer3 = wx.BoxSizer( wx.VERTICAL )

        staticBox31 = wx.StaticBoxSizer( wx.StaticBox( self.tab_domain, wx.ID_ANY, u"Windfield / Hazard Domain" ), wx.VERTICAL )

        sizer311n = wx.BoxSizer( wx.HORIZONTAL )

        self.label3111n = wx.StaticText( self.tab_domain, wx.ID_ANY, u"Country:   ", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label3111n.Wrap( -1 )
        sizer311n.Add( self.label3111n, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.choiceBoxn2 = wx.Choice(self.tab_domain, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, self.columnChoices['Countries'], 0)
        self.choiceBoxn2.SetSelection(0)
        self.choiceBoxn2.Enable(True)
        self.Bind(wx.EVT_CHOICE, self.setCountry, self.choiceBoxn2)
        sizer311n.Add( self.choiceBoxn2, 2, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.label31112n = wx.StaticText( self.tab_domain, wx.ID_ANY, u"   ", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label31112n.Wrap( -1 )
        sizer311n.Add( self.label31112n, 2, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        staticBox31.Add( sizer311n, 1, wx.EXPAND, 5 )

        
        sizer311n2 = wx.BoxSizer( wx.HORIZONTAL )

        self.label3111n2 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"Division:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label3111n2.Wrap( -1 )
        sizer311n2.Add( self.label3111n2, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.choiceBoxn2b = wx.Choice(self.tab_domain, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, [], 0)
        self.choiceBoxn2b.SetSelection(0)
        self.choiceBoxn2b.Enable(False)
        self.Bind(wx.EVT_CHOICE, self.setDivision, self.choiceBoxn2b)
        sizer311n2.Add( self.choiceBoxn2b, 4, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        staticBox31.Add( sizer311n2, 1, wx.EXPAND, 5 )


        sizer311n2b = wx.BoxSizer( wx.HORIZONTAL )

        self.label3111n2b = wx.StaticText( self.tab_domain, wx.ID_ANY, u"Location:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label3111n2b.Wrap( -1 )
        sizer311n2b.Add( self.label3111n2b, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.choiceBoxn2bb = wx.Choice(self.tab_domain, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, [], 0)
        self.choiceBoxn2bb.SetSelection(0)
        self.choiceBoxn2bb.Enable(False)
        self.Bind(wx.EVT_CHOICE, self.setLocation, self.choiceBoxn2bb)
        sizer311n2b.Add( self.choiceBoxn2bb, 4, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        staticBox31.Add( sizer311n2b, 1, wx.EXPAND, 5 )

        
        
        sizer311n3 = wx.BoxSizer( wx.HORIZONTAL )

        self.label3111n3 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"                                  Domain Limits (automatically selected):", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label3111n3.Enable(False)
        self.label3111n3.Wrap( -1 )
        sizer311n3.Add( self.label3111n3, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        staticBox31.Add( sizer311n3, 1, wx.EXPAND, 5 )
        

        sizer311 = wx.BoxSizer( wx.HORIZONTAL )

        self.label3111 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"                                  Longitude (deg E):", wx.DefaultPosition, wx.DefaultSize, 1 )
        self.label3111.Wrap( -1 )
        self.label3111.Enable(False)
        sizer311.Add( self.label3111, 2, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.label3112 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"min", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label3112.Wrap( -1 )
        self.label3112.Enable(False)
        sizer311.Add( self.label3112, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.textCtrl3113 = wx.TextCtrl( self.tab_domain, 3113, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_READONLY )
        self.textCtrl3113.Enable(False)
        sizer311.Add( self.textCtrl3113, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        self.pointers.append((self.textCtrl3113, 'Region', 'gridLimit', 'xMin'))        

        self.label3114 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"max", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label3114.Enable(False)
        self.label3114.Wrap( -1 )
        sizer311.Add( self.label3114, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.textCtrl3115 = wx.TextCtrl( self.tab_domain, 3115, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_READONLY )
        self.textCtrl3115.Enable(False)
        sizer311.Add( self.textCtrl3115, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        self.pointers.append((self.textCtrl3115, 'Region', 'gridLimit', 'xMax'))

        staticBox31.Add( sizer311, 1, wx.EXPAND, 5 )

        sizer312 = wx.BoxSizer( wx.HORIZONTAL )

        self.label3121 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"                                  Latitude (deg N):", wx.DefaultPosition, wx.DefaultSize, 1 )
        self.label3121.Wrap( -1 )
        self.label3121.Enable(False)
        sizer312.Add( self.label3121, 2, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.label3122 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"min", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label3122.Enable(False)
        self.label3122.Wrap( -1 )
        sizer312.Add( self.label3122, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.textCtrl3123 = wx.TextCtrl( self.tab_domain, 3123, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_READONLY )
        self.textCtrl3123.Enable(False)
        sizer312.Add( self.textCtrl3123, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        self.pointers.append((self.textCtrl3123, 'Region', 'gridLimit', 'yMin'))

        self.label3124 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"max", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label3124.Enable(False)
        self.label3124.Wrap( -1 )
        sizer312.Add( self.label3124, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.textCtrl3125 = wx.TextCtrl( self.tab_domain, 3125, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_READONLY )
        self.textCtrl3125.Enable(False)
        sizer312.Add( self.textCtrl3125, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        self.pointers.append((self.textCtrl3125, 'Region', 'gridLimit', 'yMax'))

        staticBox31.Add( sizer312, 1, wx.EXPAND, 5 )
        
        sizer3.Add( staticBox31, 5, wx.ALL|wx.EXPAND, 5 )

        #staticBox32 = wx.StaticBoxSizer( wx.StaticBox( self.tab_domain, wx.ID_ANY, u"Track Generator Domain" ), wx.VERTICAL )

        #sizer321 = wx.BoxSizer( wx.HORIZONTAL )

        #self.label3211 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"Longitude (deg E):", wx.DefaultPosition, wx.DefaultSize, 0 )
        #self.label3211.Wrap( -1 )
        #sizer321.Add( self.label3211, 2, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        #self.label3212 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"min", wx.DefaultPosition, wx.DefaultSize, 0 )
        #self.label3212.Wrap( -1 )
        #sizer321.Add( self.label3212, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        #self.textCtrl3213 = wx.TextCtrl( self.tab_domain, 3213, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
        #sizer321.Add( self.textCtrl3213, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        #self.pointers.append((self.textCtrl3213, 'TrackGenerator', 'gridLimit', 'xMin'))

        #self.label3214 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"max", wx.DefaultPosition, wx.DefaultSize, 0 )
        #self.label3214.Wrap( -1 )
        #sizer321.Add( self.label3214, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        #self.textCtrl3215 = wx.TextCtrl( self.tab_domain, 3215, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
        #sizer321.Add( self.textCtrl3215, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        #self.pointers.append((self.textCtrl3215, 'TrackGenerator', 'gridLimit', 'xMax'))

        #staticBox32.Add( sizer321, 1, wx.EXPAND, 5 )

        # sizer322 = wx.BoxSizer( wx.HORIZONTAL )

        # self.label3221 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"Latitude (deg N):", wx.DefaultPosition, wx.DefaultSize, 0 )
        # self.label3221.Wrap( -1 )
        # sizer322.Add( self.label3221, 2, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        # self.label3222 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"min", wx.DefaultPosition, wx.DefaultSize, 0 )
        # self.label3222.Wrap( -1 )
        # sizer322.Add( self.label3222, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        # self.textCtrl3223 = wx.TextCtrl( self.tab_domain, 3223, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
        # sizer322.Add( self.textCtrl3223, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        # self.pointers.append((self.textCtrl3223, 'TrackGenerator', 'gridLimit', 'yMin'))

        # self.label3224 = wx.StaticText( self.tab_domain, wx.ID_ANY, u"max", wx.DefaultPosition, wx.DefaultSize, 0 )
        # self.label3224.Wrap( -1 )
        # sizer322.Add( self.label3224, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        # self.textCtrl3225 = wx.TextCtrl( self.tab_domain, 3225, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
        # sizer322.Add( self.textCtrl3225, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        # self.pointers.append((self.textCtrl3225, 'TrackGenerator', 'gridLimit', 'yMax'))

        # staticBox32.Add( sizer322, 1, wx.EXPAND, 5 )

        # sizer323 = wx.BoxSizer( wx.HORIZONTAL )

        # self.button3231 = wx.Button( self.tab_domain, wx.ID_ANY, u"Autocalculate", wx.DefaultPosition, wx.DefaultSize, 0 )
        # sizer323.Add( self.button3231, 0, wx.ALL, 5 )
        # self.button3231.Bind(wx.EVT_BUTTON, self.calculateTrackGeneratorDomain)

        # self.button3232 = wx.Button( self.tab_domain, 3232, u"Help", wx.DefaultPosition, wx.DefaultSize, 0 )
        # sizer323.Add( self.button3232, 0, wx.ALL, 5 )
        # self.button3232.Bind(wx.EVT_BUTTON, self.loadHelp)

        # staticBox32.Add( sizer323, 1, wx.ALIGN_RIGHT, 5 )

        #sizer3.Add( staticBox32, 4, wx.ALL|wx.EXPAND, 5 )

        #padder33 = wx.BoxSizer( wx.VERTICAL )

        #sizer3.Add( padder33, 1, wx.EXPAND, 5 )

        self.tab_domain.SetSizer( sizer3 )

        self.tab_domain.Layout()
        sizer3.Fit( self.tab_domain )
        self.notebook_panel.AddPage( self.tab_domain, u"Domain", False )
        self.userhelptext.append(" Select a country, division and location for assessing the cyclonic wind hazard\n" + \
                                 " (e.g. Australia -> Western Australia -> Port Hedland).\n\n" + \
                                 " The domain limits will then be automatically determined.\n\n" + \
                                 " Note: locations close to the equator (i.e. within 6 degrees) or poleward of \n" + \
                                 " +/- 25 degrees latitude should not be chosen.")
        #---------------------------------------------------------


        #--------------------- Create RMW Tab --------------------
        if self.display_RMW_tab:
            self.tab_rmw = wx.Panel( self.notebook_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
            sizer9 = wx.BoxSizer( wx.VERTICAL )
            staticBox91 = wx.StaticBoxSizer( wx.StaticBox( self.tab_rmw, wx.ID_ANY, u"Radius of Maximum Winds" ), wx.VERTICAL )

            self.radioBtn911 = wx.RadioButton( self.tab_rmw, wx.ID_ANY, u"Use Log-Normal Distribution", wx.DefaultPosition, wx.DefaultSize, 0 )
            staticBox91.Add( self.radioBtn911, 2, wx.ALL, 5 )
            self.radioBtn911.Bind(wx.EVT_RADIOBUTTON, self.onRMWsel1)

            sizer912 = wx.BoxSizer( wx.HORIZONTAL )
            sizer9121 = wx.BoxSizer( wx.VERTICAL )

            sizer912.Add( sizer9121, 1, wx.EXPAND, 5 )

            self.label9122 = wx.StaticText( self.tab_rmw, wx.ID_ANY, u"mean (km)", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.label9122.Wrap( -1 )
            sizer912.Add( self.label9122, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            self.textCtrl9123 = wx.TextCtrl( self.tab_rmw, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, 0 )
            sizer912.Add( self.textCtrl9123, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.pointers.append((self.textCtrl9123, 'RMW', 'mean', 0))

            self.label9124 = wx.StaticText( self.tab_rmw, wx.ID_ANY, u"sigma", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.label9124.Wrap( -1 )
            sizer912.Add( self.label9124, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.textCtrl9125 = wx.TextCtrl( self.tab_rmw, wx.ID_ANY, u"", wx.DefaultPosition, wx.DefaultSize, wx.TE_PROCESS_ENTER )
            sizer912.Add( self.textCtrl9125, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.pointers.append((self.textCtrl9125, 'RMW', 'sigma', 0))

            staticBox91.Add( sizer912, 2, wx.EXPAND, 5 )

            self.radioBtn913 = wx.RadioButton( self.tab_rmw, wx.ID_ANY, u"Estimate Distribution from Input Data  (Note: this requires RMW column)", wx.DefaultPosition, wx.DefaultSize, 0 )
            staticBox91.Add( self.radioBtn913, 2, wx.ALL, 5 )
            self.radioBtn913.Bind(wx.EVT_RADIOBUTTON, self.onRMWsel2)

            self.m_button914 = wx.Button( self.tab_rmw, 914, u"Help", wx.DefaultPosition, wx.DefaultSize, 0 )
            staticBox91.Add( self.m_button914, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT|wx.ALL, 5 )
            self.m_button914.Bind(wx.EVT_BUTTON, self.loadHelp)

            sizer9.Add( staticBox91, 3, wx.ALL|wx.EXPAND, 5 )
            padder92 = wx.BoxSizer( wx.VERTICAL )
            sizer9.Add( padder92, 2, wx.EXPAND, 5 )

            self.tab_rmw.SetSizer( sizer9 )
            self.tab_rmw.Layout()
            sizer9.Fit( self.tab_rmw )
            self.notebook_panel.AddPage( self.tab_rmw, u"RMW", False )
            self.userhelptext.append(" Choose a Radius of Maximum Wind (RMW) distribution.  \n Two options are available: \n (1) lognormal distribution \n (2) distribution determined from the input data using kernel density estimation.")
        #---------------------------------------------------------


        #-------------------- Create MSLP Tab --------------------
        if self.display_MSLP_tab:
            self.tab_mslp = wx.Panel( self.notebook_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
            sizer4 = wx.BoxSizer( wx.VERTICAL )

            staticBox41 = wx.StaticBoxSizer( wx.StaticBox( self.tab_mslp, wx.ID_ANY, u"Mean Sea Level Pressure" ), wx.VERTICAL )

            self.radioBtn411 = wx.RadioButton( self.tab_mslp, wx.ID_ANY, u"Use Yearly Average", wx.DefaultPosition, wx.DefaultSize, 0 )
            staticBox41.Add( self.radioBtn411, 2, wx.ALL, 5 )
            self.radioBtn411.Bind(wx.EVT_RADIOBUTTON, self.OnYearlyAverage)

            self.radioBtn412 = wx.RadioButton( self.tab_mslp, wx.ID_ANY, u"Use Seasonal Average", wx.DefaultPosition, wx.DefaultSize, 0 )
            staticBox41.Add( self.radioBtn412, 2, wx.ALL, 5 )
            self.radioBtn412.Bind(wx.EVT_RADIOBUTTON, self.OnSeasonalAverage)

            sizer413 = wx.BoxSizer( wx.HORIZONTAL )

            sizer4131 = wx.BoxSizer( wx.VERTICAL )

            sizer413.Add( sizer4131, 1, wx.EXPAND, 5 )

            sizer4132 = wx.BoxSizer( wx.VERTICAL )

            self.label41321 = wx.StaticText( self.tab_mslp, wx.ID_ANY, u"Months for Averaging:", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.label41321.Wrap( -1 )
            self.label41321.Enable( False )

            sizer4132.Add( self.label41321, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            gSizer41322 = wx.GridSizer( 2, 6, 0, 0 )

            # Create month selection boxes for MSLP seasonal average
            month_strs = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
            self.month_sel = []
            for k in range(12):
                self.month_sel.append(wx.CheckBox(self.tab_mslp, wx.ID_ANY, month_strs[k], wx.DefaultPosition, wx.DefaultSize, 0))
                self.month_sel[k].Bind(wx.EVT_CHECKBOX, self.onMSLPcheckbox)
                gSizer41322.Add( self.month_sel[k], 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            sizer4132.Add( gSizer41322, 2, wx.EXPAND, 5 )

            sizer413.Add( sizer4132, 3, wx.EXPAND, 5 )

            staticBox41.Add( sizer413, 5, wx.EXPAND, 5 )

            self.button414 = wx.Button( self.tab_mslp, 414, u"Help", wx.DefaultPosition, wx.DefaultSize, 0 )
            staticBox41.Add( self.button414, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT|wx.ALL, 5 )
            self.button414.Bind(wx.EVT_BUTTON, self.loadHelp)

            sizer4.Add( staticBox41, 4, wx.ALL|wx.EXPAND, 5 )

            padder42 = wx.BoxSizer( wx.VERTICAL )

            sizer4.Add( padder42, 1, wx.EXPAND, 5 )

            self.tab_mslp.SetSizer( sizer4 )
            self.tab_mslp.Layout()
            sizer4.Fit( self.tab_mslp )
            self.notebook_panel.AddPage( self.tab_mslp, u"MSLP", False )
            self.userhelptext.append(" Choose a yearly or seasonal average Mean Sea Level Pressure climatology field.\n\n For more information, please press the 'Help' button.")
        #---------------------------------------------------------


        #--------------- Create Stats Interface Tab --------------
        if self.display_stats_tab:
            self.tab_stats = wx.Panel( self.notebook_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
            sizer5 = wx.BoxSizer( wx.VERTICAL )

            staticBox51 = wx.StaticBoxSizer( wx.StaticBox( self.tab_stats, wx.ID_ANY, u"Distribution Types" ), wx.VERTICAL )

            sizer511 = wx.BoxSizer( wx.HORIZONTAL )

            self.label5111 = wx.StaticText( self.tab_stats, wx.ID_ANY, u"1D KDE Type:", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.label5111.Wrap( -1 )
            sizer511.Add( self.label5111, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            self.radioBtn5112 = wx.RadioButton( self.tab_stats, wx.ID_ANY, self.KDETypes[0], wx.DefaultPosition, wx.DefaultSize, wx.RB_GROUP)
            sizer511.Add( self.radioBtn5112, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.radioBtn5112.Bind(wx.EVT_RADIOBUTTON, self.updateKDE1DSettings)

            self.radioBtn5113 = wx.RadioButton( self.tab_stats, wx.ID_ANY, self.KDETypes[1], wx.DefaultPosition, wx.DefaultSize, 0 )
            sizer511.Add( self.radioBtn5113, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.radioBtn5113.Bind(wx.EVT_RADIOBUTTON, self.updateKDE1DSettings)

            self.radioBtn5114 = wx.RadioButton( self.tab_stats, wx.ID_ANY, self.KDETypes[2], wx.DefaultPosition, wx.DefaultSize, 0 )
            sizer511.Add( self.radioBtn5114, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.radioBtn5114.Bind(wx.EVT_RADIOBUTTON, self.updateKDE1DSettings)

            self.radioBtn5115 = wx.RadioButton( self.tab_stats, wx.ID_ANY, self.KDETypes[3], wx.DefaultPosition, wx.DefaultSize, 0 )
            sizer511.Add( self.radioBtn5115, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.radioBtn5115.Bind(wx.EVT_RADIOBUTTON, self.updateKDE1DSettings)

            staticBox51.Add( sizer511, 1, wx.EXPAND, 5 )

            sizer512 = wx.BoxSizer( wx.HORIZONTAL )

            self.label5121 = wx.StaticText( self.tab_stats, wx.ID_ANY, u"2D KDE Type:", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.label5121.Wrap( -1 )
            sizer512.Add( self.label5121, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            self.radioBtn5122 = wx.RadioButton( self.tab_stats, wx.ID_ANY, self.KDETypes[0], wx.DefaultPosition, wx.DefaultSize, wx.RB_GROUP)
            sizer512.Add( self.radioBtn5122, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.radioBtn5122.Bind(wx.EVT_RADIOBUTTON, self.updateKDE2DSettings)

            self.radioBtn5123 = wx.RadioButton( self.tab_stats, wx.ID_ANY, self.KDETypes[1], wx.DefaultPosition, wx.DefaultSize, 0 )
            sizer512.Add( self.radioBtn5123, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.radioBtn5123.Bind(wx.EVT_RADIOBUTTON, self.updateKDE2DSettings)

            self.radioBtn5124 = wx.RadioButton( self.tab_stats, wx.ID_ANY, self.KDETypes[2], wx.DefaultPosition, wx.DefaultSize, 0 )
            sizer512.Add( self.radioBtn5124, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.radioBtn5124.Bind(wx.EVT_RADIOBUTTON, self.updateKDE2DSettings)

            self.radioBtn5125 = wx.RadioButton( self.tab_stats, wx.ID_ANY, self.KDETypes[3], wx.DefaultPosition, wx.DefaultSize, 0 )
            sizer512.Add( self.radioBtn5125, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.radioBtn5125.Bind(wx.EVT_RADIOBUTTON, self.updateKDE2DSettings)

            staticBox51.Add( sizer512, 1, wx.EXPAND, 5 )

            sizer5.Add( staticBox51, 4, wx.ALL|wx.EXPAND, 5 )

            staticBox52 = wx.StaticBoxSizer( wx.StaticBox( self.tab_stats, wx.ID_ANY, u"Box Size for Distribution Fitting" ), wx.VERTICAL )

            sizer521 = wx.BoxSizer( wx.HORIZONTAL )

            self.label5211 = wx.StaticText( self.tab_stats, wx.ID_ANY, u"Grid Space:", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.label5211.Wrap( -1 )
            sizer521.Add( self.label5211, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            self.label5212 = wx.StaticText( self.tab_stats, wx.ID_ANY, u"deg lon", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.label5212.Wrap( -1 )
            sizer521.Add( self.label5212, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            self.textCtrl5213 = wx.TextCtrl( self.tab_stats, 5213, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
            sizer521.Add( self.textCtrl5213, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.pointers.append((self.textCtrl5213, 'StatInterface', 'gridSpace', 'x'))

            self.label5214 = wx.StaticText( self.tab_stats, wx.ID_ANY, u"deg lat", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.label5214.Wrap( -1 )
            sizer521.Add( self.label5214, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            self.textCtrl5215 = wx.TextCtrl( self.tab_stats, 5215, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
            sizer521.Add( self.textCtrl5215, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.pointers.append((self.textCtrl5215, 'StatInterface', 'gridSpace', 'y'))

            staticBox52.Add( sizer521, 1, wx.EXPAND, 5 )

            sizer522 = wx.BoxSizer( wx.HORIZONTAL )

            self.label5221 = wx.StaticText( self.tab_stats, wx.ID_ANY, u"Grid Expansion Increment:", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.label5221.Wrap( -1 )
            sizer522.Add( self.label5221, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            self.label5222 = wx.StaticText( self.tab_stats, wx.ID_ANY, u"deg lon", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.label5222.Wrap( -1 )
            sizer522.Add( self.label5222, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            self.textCtrl5223 = wx.TextCtrl( self.tab_stats, 5223, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
            sizer522.Add( self.textCtrl5223, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.pointers.append((self.textCtrl5223, 'StatInterface', 'gridInc', 'x'))

            self.label5224 = wx.StaticText( self.tab_stats, wx.ID_ANY, u"deg lat", wx.DefaultPosition, wx.DefaultSize, 0)
            self.label5224.Wrap( -1 )
            sizer522.Add( self.label5224, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            self.textCtrl5225 = wx.TextCtrl( self.tab_stats, 5225, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_PROCESS_ENTER)
            sizer522.Add( self.textCtrl5225, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.pointers.append((self.textCtrl5225, 'StatInterface', 'gridInc', 'y'))

            staticBox52.Add( sizer522, 1, wx.EXPAND, 5 )

            sizer5.Add( staticBox52, 5, wx.ALL|wx.EXPAND, 5 )

            padder53 = wx.BoxSizer( wx.VERTICAL )

            sizer5.Add( padder53, 2, wx.EXPAND, 5 )

            self.tab_stats.SetSizer( sizer5 )
            self.tab_stats.Layout()
            sizer5.Fit( self.tab_stats )
            self.notebook_panel.AddPage( self.tab_stats, u"Stats Interface", False )
            self.userhelptext.append(" Storm statistics are determined at each location by applying kernel density \n estimation.  The settings governing this process can be modified above.\n" + \
                                     " The 'Grid Space' setting determines the size of the grid used for aggregating \n samples.  If insufficient samples are found, then the grid will expand in increments \n specified by the 'Grid Expansion Increment' setting.")
        #---------------------------------------------------------


        #--------------- Create Track Generator Tab --------------
        self.tab_trackgen = wx.Panel( self.notebook_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.tab_trackgen.SetMinSize( wx.Size( -1,1 ) )

        sizer6 = wx.BoxSizer( wx.VERTICAL )

        staticBox61 = wx.StaticBoxSizer( wx.StaticBox( self.tab_trackgen, wx.ID_ANY, u"Simulation Settings" ), wx.VERTICAL )

        sizer611 = wx.BoxSizer( wx.HORIZONTAL )

        self.label6111 = wx.StaticText( self.tab_trackgen, wx.ID_ANY, u"Number of Years to Simulate:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label6111.Wrap( -1 )
        sizer611.Add( self.label6111, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.textCtrl6112 = wx.TextCtrl( self.tab_trackgen, 6112, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
        sizer611.Add( self.textCtrl6112, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        self.pointers.append((self.textCtrl6112, 'TrackGenerator', 'NumSimulations', 0))

        staticBox61.Add( sizer611, 1, wx.EXPAND, 5 )

        #sizer612 = wx.BoxSizer( wx.HORIZONTAL )

        #self.label6121 = wx.StaticText( self.tab_trackgen, wx.ID_ANY, u"Years per Simulation:", wx.DefaultPosition, wx.DefaultSize, 0 )
        #self.label6121.Wrap( -1 )
        #sizer612.Add( self.label6121, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        #self.textCtrl6122 = wx.TextCtrl( self.tab_trackgen, 6122, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
        #sizer612.Add( self.textCtrl6122, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        #self.pointers.append((self.textCtrl6122, 'TrackGenerator', 'YearsPerSimulation', 0))

        #staticBox61.Add( sizer612, 1, wx.EXPAND, 5 )

        sizer6.Add( staticBox61, 2, wx.ALL|wx.EXPAND, 5 )

        #staticBox62 = wx.StaticBoxSizer( wx.StaticBox( self.tab_trackgen, wx.ID_ANY, u"Storm Frequency" ), wx.VERTICAL )

        #sizer621 = wx.BoxSizer( wx.HORIZONTAL )

        #self.label6211 = wx.StaticText( self.tab_trackgen, wx.ID_ANY, u"Annual Frequency in Track Generator Domain:", wx.DefaultPosition, wx.DefaultSize, 0 )
        #self.label6211.Wrap( -1 )
        #sizer621.Add( self.label6211, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        #self.textCtrl6212 = wx.TextCtrl( self.tab_trackgen, 6212, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
        #sizer621.Add( self.textCtrl6212, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        #self.pointers.append((self.textCtrl6212, 'TrackGenerator', 'Frequency', 0))

        #staticBox62.Add( sizer621, 1, wx.EXPAND, 5 )

        #sizer622 = wx.BoxSizer( wx.HORIZONTAL )

        #self.button6221 = wx.Button( self.tab_trackgen, wx.ID_ANY, u"Autocalculate", wx.DefaultPosition, wx.DefaultSize, 0 )
        #sizer622.Add( self.button6221, 0, wx.ALL, 5 )
        #self.button6221.Bind(wx.EVT_BUTTON, self.calculateFrequency)

        #self.button6222 = wx.Button( self.tab_trackgen, 6222, u"Help", wx.DefaultPosition, wx.DefaultSize, 0 )
        #sizer622.Add( self.button6222, 0, wx.ALL, 5 )
        #self.button6222.Bind(wx.EVT_BUTTON, self.loadHelp)

        #staticBox62.Add( sizer622, 1, wx.ALIGN_RIGHT, 5 )

        #sizer6.Add( staticBox62, 2, wx.ALL|wx.EXPAND, 5 )

        padder63 = wx.BoxSizer( wx.VERTICAL )

        sizer6.Add( padder63, 4, wx.EXPAND, 5 )

        self.tab_trackgen.SetSizer( sizer6 )
        self.tab_trackgen.Layout()
        sizer6.Fit( self.tab_trackgen )
        self.notebook_panel.AddPage( self.tab_trackgen, u"Track Generator", False )
        self.userhelptext.append(" Choose the number of years to simulate.\n\n" + \
                                 " This should be much larger than the return period of interest \n" + \
                                 " (e.g. 5000 year simulation for a 500 year return period).")
        #---------------------------------------------------------


        #------------------ Create Windfield Tab -----------------
        if self.display_windfield_tab:
            self.tab_windfield = wx.Panel( self.notebook_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
            sizer7 = wx.BoxSizer( wx.VERTICAL )

            staticBox71 = wx.StaticBoxSizer( wx.StaticBox( self.tab_windfield, wx.ID_ANY, u"WindField Profile" ), wx.VERTICAL )

            sizer711 = wx.BoxSizer( wx.HORIZONTAL )

            self.radioBtn7111 = wx.RadioButton( self.tab_windfield, wx.ID_ANY, self.profileTypes_display[0], wx.DefaultPosition, wx.DefaultSize, wx.RB_GROUP)
            sizer711.Add( self.radioBtn7111, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.radioBtn7111.Bind(wx.EVT_RADIOBUTTON, self.onOtherProfile)

            self.radioBtn7112 = wx.RadioButton( self.tab_windfield, wx.ID_ANY, self.profileTypes_display[1], wx.DefaultPosition, wx.DefaultSize, 0 )
            sizer711.Add( self.radioBtn7112, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.radioBtn7112.Bind(wx.EVT_RADIOBUTTON, self.onHollandProfile)

            # self.radioBtn7113 = wx.RadioButton( self.tab_windfield, wx.ID_ANY, self.profileTypes_display[2], wx.DefaultPosition, wx.DefaultSize, 0 )
            # sizer711.Add( self.radioBtn7113, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            # self.radioBtn7113.Bind(wx.EVT_RADIOBUTTON, self.onHollandProfile)

            # self.radioBtn7114 = wx.RadioButton( self.tab_windfield, wx.ID_ANY, self.profileTypes_display[3], wx.DefaultPosition, wx.DefaultSize, 0 )
            # sizer711.Add( self.radioBtn7114, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            # self.radioBtn7114.Bind(wx.EVT_RADIOBUTTON, self.onOtherProfile)

            #sizer7119 = wx.BoxSizer( wx.HORIZONTAL )

            # self.radioBtn7115 = wx.RadioButton( self.tab_windfield, wx.ID_ANY, self.profileTypes_display[4], wx.DefaultPosition, wx.DefaultSize, 0 )
            # sizer7119.Add( self.radioBtn7115, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            # self.radioBtn7115.Bind(wx.EVT_RADIOBUTTON, self.onOtherProfile)

            # self.radioBtn7116 = wx.RadioButton( self.tab_windfield, wx.ID_ANY, self.profileTypes_display[5], wx.DefaultPosition, wx.DefaultSize, 0 )
            # sizer7119.Add( self.radioBtn7116, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            # self.radioBtn7116.Bind(wx.EVT_RADIOBUTTON, self.onOtherProfile)

            staticBox71.Add( sizer711, 1, wx.EXPAND, 5 )
            #staticBox71.Add( sizer7119, 1, wx.EXPAND, 5 )

            sizer712 = wx.BoxSizer( wx.HORIZONTAL )

            self.label7121 = wx.StaticText( self.tab_windfield, wx.ID_ANY, u"Beta:", wx.DefaultPosition, wx.DefaultSize, 0 )
            self.label7121.Wrap( -1 )
            sizer712.Add( self.label7121, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            self.textCtrl7122 = wx.TextCtrl( self.tab_windfield, 7122, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_PROCESS_ENTER )
            sizer712.Add( self.textCtrl7122, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            self.pointers.append((self.textCtrl7122, 'WindfieldInterface', 'beta', 0))

            #self.label7123 = wx.StaticText( self.tab_windfield, wx.ID_ANY, u"Beta1:", wx.DefaultPosition, wx.DefaultSize, 0 )
            #self.label7123.Wrap( -1 )

            #sizer712.Add( self.label7123, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            #self.textCtrl7124 = wx.TextCtrl( self.tab_windfield, 7124, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
            #sizer712.Add( self.textCtrl7124, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            #self.pointers.append((self.textCtrl7124, 'WindfieldInterface', 'beta1', 0))

            #self.label7125 = wx.StaticText( self.tab_windfield, wx.ID_ANY, u"Beta2:", wx.DefaultPosition, wx.DefaultSize, 0 )
            #self.label7125.Wrap( -1 )

            #sizer712.Add( self.label7125, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            #self.textCtrl7126 = wx.TextCtrl( self.tab_windfield, 7126, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_PROCESS_ENTER )
            #sizer712.Add( self.textCtrl7126, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            #self.pointers.append((self.textCtrl7126, 'WindfieldInterface', 'beta2', 0))

            staticBox71.Add( sizer712, 1, wx.EXPAND, 5 )

            sizer7.Add( staticBox71, 2, wx.ALL|wx.EXPAND, 5 )

            #staticBox72 = wx.StaticBoxSizer( wx.StaticBox( self.tab_windfield, wx.ID_ANY, u"Boundary Layer Model" ), wx.HORIZONTAL )

            #self.radioBtn721 = wx.RadioButton( self.tab_windfield, wx.ID_ANY, self.windFieldTypes_display[0], wx.DefaultPosition, wx.DefaultSize, wx.RB_GROUP)
            #staticBox72.Add( self.radioBtn721, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            #self.radioBtn721.Bind(wx.EVT_RADIOBUTTON, self.onKepert)

            #self.radioBtn722 = wx.RadioButton( self.tab_windfield, wx.ID_ANY, self.windFieldTypes_display[1], wx.DefaultPosition, wx.DefaultSize, 0 )
            #staticBox72.Add( self.radioBtn722, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            #self.radioBtn722.Bind(wx.EVT_RADIOBUTTON, self.onMcConochieOrHubbert)

            #self.radioBtn723 = wx.RadioButton( self.tab_windfield, wx.ID_ANY, self.windFieldTypes_display[2], wx.DefaultPosition, wx.DefaultSize, 0 )
            #staticBox72.Add( self.radioBtn723, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            #self.radioBtn723.Bind(wx.EVT_RADIOBUTTON, self.onMcConochieOrHubbert)

            #self.label724 = wx.StaticText( self.tab_windfield, wx.ID_ANY, u"thetaMax:", wx.DefaultPosition, wx.DefaultSize, 0 )
            #self.label724.Wrap( -1 )
            #staticBox72.Add( self.label724, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            #self.textCtrl725 = wx.TextCtrl( self.tab_windfield, 725, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
            #staticBox72.Add( self.textCtrl725, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            #self.pointers.append((self.textCtrl725, 'WindfieldInterface', 'thetaMax', 0))

            #sizer7.Add( staticBox72, 2, wx.ALL|wx.EXPAND, 5 )

            #staticBox73 = wx.StaticBoxSizer( wx.StaticBox( self.tab_windfield, wx.ID_ANY, u"Resolution" ), wx.HORIZONTAL )

            #self.label731 = wx.StaticText( self.tab_windfield, wx.ID_ANY, u"Windfield resolution (in degrees):", wx.DefaultPosition, wx.DefaultSize, 0 )
            #self.label731.Wrap( -1 )
            #staticBox73.Add( self.label731, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

            #self.textCtrl732 = wx.TextCtrl( self.tab_windfield, 732, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_PROCESS_ENTER )
            #staticBox73.Add( self.textCtrl732, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
            #self.pointers.append((self.textCtrl732, 'WindfieldInterface', 'Resolution', 0))

            #sizer7.Add( staticBox73, 2, wx.ALL|wx.EXPAND, 5 )
            padder74 = wx.BoxSizer( wx.VERTICAL )
            sizer7.Add( padder74, 3, wx.EXPAND, 5 )

            self.tab_windfield.SetSizer( sizer7 )
            self.tab_windfield.Layout()
            sizer7.Fit( self.tab_windfield )
            self.notebook_panel.AddPage( self.tab_windfield, u"Windfield", False )
            self.userhelptext.append(" Select a windfield profile.")
        #---------------------------------------------------------


        #------------------- Create Hazard Tab -------------------
        self.tab_hazard = wx.Panel( self.notebook_panel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        sizer8 = wx.BoxSizer( wx.VERTICAL )

        staticBox81 = wx.StaticBoxSizer( wx.StaticBox( self.tab_hazard, wx.ID_ANY, u"Return Periods" ), wx.HORIZONTAL )

        self.label811 = wx.StaticText( self.tab_hazard, wx.ID_ANY, u"Years:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label811.Wrap( -1 )
        staticBox81.Add( self.label811, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        self.textCtrl812 = wx.TextCtrl( self.tab_hazard, 812, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_PROCESS_ENTER )
        staticBox81.Add( self.textCtrl812, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        self.pointers.append((self.textCtrl812, 'HazardInterface', 'Years', 0))

        sizer8.Add( staticBox81, 1, wx.ALL|wx.EXPAND, 5 )

        staticBox82 = wx.StaticBoxSizer( wx.StaticBox( self.tab_hazard, wx.ID_ANY, u"Plotting Units" ), wx.HORIZONTAL )
        self.label8221 = wx.StaticText( self.tab_hazard, wx.ID_ANY, u"Speed Units:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label8221.Wrap( -1 )
        staticBox82.Add( self.label8221, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 6 )

        self.choiceBoxn822 = wx.Choice(self.tab_hazard, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, self.speedUnitsLongname, 0)
        self.choiceBoxn822.SetSelection(0)
        self.Bind(wx.EVT_CHOICE, self.setPlottingUnits, self.choiceBoxn822)
        staticBox82.Add( self.choiceBoxn822, 2, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 6 )

        self.label823 = wx.StaticText( self.tab_hazard, wx.ID_ANY, u"  ", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.label823.Wrap( -1 )
        staticBox82.Add( self.label823, 3, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 6 )
        
        sizer8.Add( staticBox82, 1, wx.ALL|wx.EXPAND, 5 )

        padder82 = wx.BoxSizer( wx.VERTICAL )
        sizer8.Add( padder82, 1, wx.EXPAND, 5 )        

        self.tab_hazard.SetSizer( sizer8 )
        self.tab_hazard.Layout()
        sizer8.Fit( self.tab_hazard )
        self.notebook_panel.AddPage( self.tab_hazard, u"Hazard", False )
        self.userhelptext.append(" Enter the required return periods for the wind hazard calculation. \n The input should be a list of years delimited by commas.\n\n Also select the units for plotting the hazard results.")
        #---------------------------------------------------------


        staticBox02.Add( self.notebook_panel, 1, wx.ALL|wx.EXPAND, 5 )

        sizer0.Add( staticBox02, 35, wx.ALL|wx.EXPAND, 5 )

        sbSizer16 = wx.StaticBoxSizer( wx.StaticBox( self.backpanel, wx.ID_ANY, u"User Hints and Info" ), wx.HORIZONTAL )

        self.m_scrolledWindow2 = wx.ScrolledWindow(self.backpanel, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.VSCROLL|wx.SUNKEN_BORDER)
        self.m_scrolledWindow2.SetBackgroundColour(wx.Colour(255, 255, 255));
        self.m_scrolledWindow2.SetScrollRate(5, 5)
        self.bSizer421 = wx.BoxSizer(wx.VERTICAL)

        self.m_staticText36 = wx.StaticText( self.m_scrolledWindow2, wx.ID_ANY, initialmsg, wx.DefaultPosition, wx.DefaultSize, wx.EXPAND )
        self.m_staticText36.Wrap( -1 )

        self.bSizer421.Add( self.m_staticText36, 1, wx.ALL|wx.EXPAND, 0 )

        self.m_scrolledWindow2.SetSizer( self.bSizer421 )
        self.m_scrolledWindow2.Layout()
        self.bSizer421.Fit( self.m_scrolledWindow2 )
        sbSizer16.Add( self.m_scrolledWindow2, 1, wx.EXPAND |wx.ALL, 5 )

        bSizer43 = wx.BoxSizer( wx.VERTICAL )

        self.m_button8 = wx.Button( self.backpanel, wx.ID_ANY, u"Save", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer43.Add( self.m_button8, 0, wx.ALIGN_BOTTOM|wx.ALL|wx.EXPAND, 5 )
        self.m_button8.Bind(wx.EVT_BUTTON, self.menu_sel_file_save)
        self.m_button8.Enable(False)

        self.m_button9 = wx.Button( self.backpanel, wx.ID_ANY, u"Save && Continue", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer43.Add( self.m_button9, 0, wx.ALIGN_BOTTOM|wx.ALL|wx.EXPAND, 5 )
        self.m_button9.Bind(wx.EVT_BUTTON, self.save_and_continue)
        self.m_button9.Enable(False)

        sbSizer16.Add( bSizer43, 0, wx.ALIGN_CENTER, 5 )
        sizer0.Add( sbSizer16, 10, wx.ALL|wx.EXPAND, 5 )

        self.backpanel.SetSizer( sizer0 )
        self.backpanel.Layout()
        sizer0.Fit( self.backpanel )
        sizer00.Add( self.backpanel, 1, wx.EXPAND|wx.TOP, 1 )

        self.SetSizer( sizer00 )

        self.StatusBar = self.CreateStatusBar( 1, wx.ST_SIZEGRIP, wx.ID_ANY )
        self.scrolledWindow122.Enable(True)

        self.idMap = {}
        self.initialiseSettings()
        self.refreshGUIsettings()

        for k in range(len(self.pointers)):
            self.pointers[k][0].Bind(wx.EVT_TEXT_ENTER, self.checkValue)
            self.pointers[k][0].Bind(wx.EVT_KILL_FOCUS, self.checkValue)
            self.idMap[self.pointers[k][0].Id] = k

        self.notebook_panel.Enable(False)
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.onTabChange)
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.onTabChanged)
        self.Layout()


    def initialiseSettings(self):
        self.configDict = {}
        self.configDictSections = []
        self.nextRecNo = 0
        emptyDomain = "{'xMin': ,'xMax': ,'yMin': ,'yMax': }"
        self.setConfigKeyValue('Actions', '; TCRM modules to execute', '', valueDataType='str')
        self.setConfigKeyValue('Actions', 'DataProcess', 'True', valueDataType='logical', allowOverwrite=False)
        self.setConfigKeyValue('Actions', 'ExecuteStat', 'True', valueDataType='logical', allowOverwrite=False)
        self.setConfigKeyValue('Actions', 'ExecuteTrackGenerator', 'True', valueDataType='logical', allowOverwrite=False)
        self.setConfigKeyValue('Actions', 'ExecuteWindfield', 'True', valueDataType='logical', allowOverwrite=False)
        self.setConfigKeyValue('Actions', 'ExecuteHazard', 'True', valueDataType='logical', allowOverwrite=False)
        self.setConfigKeyValue('Actions', 'PlotHazard', 'True', valueDataType='logical', allowOverwrite=False)
        self.setConfigKeyValue('DataProcess', 'InputFile', os.path.join(self.tcrmPath, 'input' + os.sep), valueDataType='str')
        self.setConfigKeyValue('DataProcess', 'Source', self.defaultSourceName, valueDataType='str')
        self.setConfigKeyValue('DataProcess', 'StartSeason', '1981', valueDataType='int')
        self.setConfigKeyValue('Region', '; Domain for windfield and hazard calculation', '', valueDataType='str')
        self.setConfigKeyValue('Region', 'gridLimit', str(emptyDomain), valueDataType='domain')
        self.setConfigKeyValue('Region', 'LocalityID', '-99999', valueDataType='int')
        #self.setConfigKeyValue('StatInterface', 'gridLimit', str(emptyDomain), valueDataType='domain')
        self.setConfigKeyValue('StatInterface', 'kdeType', self.KDETypes[0], valueDataType='choicelist', choicelist=self.KDETypes)
        self.setConfigKeyValue('StatInterface', 'kde2DType', self.KDETypes[3], valueDataType='choicelist', choicelist=self.KDETypes)
        self.setConfigKeyValue('StatInterface', 'kdeStep', '0.2', valueDataType='float', minVal=0)
        self.setConfigKeyValue('StatInterface', 'gridSpace', "{'x':1, 'y':1}", valueDataType='xystep')
        self.setConfigKeyValue('StatInterface', 'gridInc', "{'x':1, 'y':0.5}", valueDataType='xystep')
        #self.setConfigKeyValue('TrackGenerator', 'gridLimit', str(emptyDomain), valueDataType='domain')
        self.setConfigKeyValue('TrackGenerator', 'NumSimulations', '5000', valueDataType='int', minVal=1000, maxVal=9999)
        self.setConfigKeyValue('TrackGenerator', 'YearsPerSimulation', '1', valueDataType='int', minVal=1, maxVal=100)
        self.setConfigKeyValue('TrackGenerator', 'Frequency', '', valueDataType='float', minVal=0)
        self.setConfigKeyValue('TrackGenerator', 'gridSpace', "{'x': , 'y': }", valueDataType='xystep')
        self.setConfigKeyValue('TrackGenerator', 'gridInc', "{'x': , 'y': }", valueDataType='xystep')
        self.setConfigKeyValue('WindfieldInterface', 'NumberofFiles', '', valueDataType='int', minVal=1, maxVal=9999)
        self.setConfigKeyValue('WindfieldInterface', 'TrackPath', '{output directory}/tracks', valueDataType='str', allowOverwrite=False)
        self.setConfigKeyValue('WindfieldInterface', 'Margin', '2', valueDataType='float', minVal=0)
        self.setConfigKeyValue('WindfieldInterface', 'Resolution', '0.05', valueDataType='float', minVal=0)
        self.setConfigKeyValue('WindfieldInterface', 'Source', 'TCRM', valueDataType='str')
        self.setConfigKeyValue('WindfieldInterface', 'profileType', self.profileTypes[0], valueDataType='choicelist', choicelist=self.profileTypes)
        self.setConfigKeyValue('WindfieldInterface', 'windFieldType', self.windFieldTypes[0], valueDataType='choicelist', choicelist=self.windFieldTypes)
        self.setConfigKeyValue('WindfieldInterface', 'beta', '1.5', valueDataType='float', minVal=0)
        self.setConfigKeyValue('WindfieldInterface', 'beta1', '1.5', valueDataType='float', minVal=0)
        self.setConfigKeyValue('WindfieldInterface', 'beta2', '1.4', valueDataType='float', minVal=0)
        self.setConfigKeyValue('WindfieldInterface', 'thetaMax', '70.0', valueDataType='float', minVal=-180, maxVal=180)
        self.setConfigKeyValue('HazardInterface', '; Years to calculate return period wind speeds', '', valueDataType='str')
        self.setConfigKeyValue('HazardInterface', 'InputPath', '{output directory}/windfield', valueDataType='str', allowOverwrite=False)
        self.setConfigKeyValue('HazardInterface', 'Resolution', '', valueDataType='float', minVal=0)
        self.setConfigKeyValue('HazardInterface', 'Years', '25,50,100,200,250,500,1000,2000,2500', valueDataType='intlist', minVal=1, maxVal=9999)
        self.setConfigKeyValue('HazardInterface', 'NumSim', '', valueDataType='int', minVal=1, maxVal=9999)
        self.setConfigKeyValue('HazardInterface', 'MinimumRecords', '50', valueDataType='int', minVal=1, maxVal=999)
        self.setConfigKeyValue('HazardInterface', 'CalculateCI', 'False', valueDataType='logical')
        self.setConfigKeyValue('HazardInterface', 'PlotSpeedUnits', 'mps', valueDataType='choicelist', choicelist=self.speedUnits)
        self.setConfigKeyValue('Input', 'MSLPGrid', '[1 2 3 4 5 6 7 8 9 10 11 12]', valueDataType='intlist', minVal=1, maxVal=12)
        self.setConfigKeyValue('Input', 'LandMask', '{tcrm directory}/input/landmask.nc', valueDataType='str', allowOverwrite=False)
        self.setConfigKeyValue('Output', 'Path', os.path.join(self.tcrmPath, 'output'), valueDataType='str')
        self.setConfigKeyValue('Logging', 'LogFile', '{output directory}/log/{filename}.log', valueDataType='str', allowOverwrite=False)
        self.setConfigKeyValue('Logging', 'LogLevel', 'INFO', valueDataType='str')
        self.setConfigKeyValue('Logging', 'Verbose', 'False', valueDataType='logical')
        self.setConfigKeyValue('Process', 'DatFile', '{output directory}/process/dat/{filename}.dat', valueDataType='str', allowOverwrite=False)
        self.setConfigKeyValue('Process', 'ExcludePastProcessed', 'True', valueDataType='logical')
        self.setConfigKeyValue('RMW', 'GetRMWDistFromInputData', 'False', valueDataType='logical')
        self.setConfigKeyValue('RMW', 'mean', '50.0', valueDataType='float', minVal=0)
        self.setConfigKeyValue('RMW', 'sigma', '0.6', valueDataType='float', minVal=0)
        self.setConfigKeyValue('TCRM', '; Output track files settings', '', valueDataType='str')
        self.setConfigKeyValue('TCRM', 'Columns', 'index,age,lon,lat,speed,bearing,pressure,penv,rmax', valueDataType='str')
        self.setConfigKeyValue('TCRM', 'FieldDelimiter', ',', valueDataType='str')
        self.setConfigKeyValue('TCRM', 'NumberOfHeadingLines', '1', valueDataType='int')
        self.setConfigKeyValue('TCRM', 'SpeedUnits', 'kph', valueDataType='str')
        self.setConfigKeyValue('TCRM', 'PressureUnits', 'hPa', valueDataType='str')
        self.setConfigKeyValue(self.defaultSourceName, '; Input data file settings', '', valueDataType='str')
        self.setConfigKeyValue(self.defaultSourceName, 'Columns', '', valueDataType='str')
        self.setConfigKeyValue(self.defaultSourceName, 'FieldDelimiter', ',', valueDataType='str')
        self.setConfigKeyValue(self.defaultSourceName, 'NumberOfHeadingLines', '0', valueDataType='int')
        #self.setConfigKeyValue(self.defaultSourceName, 'SpeedUnits', 'mps', valueDataType='choicelist', choicelist=self.speedUnits)
        self.setConfigKeyValue(self.defaultSourceName, 'PressureUnits', 'hPa', valueDataType='choicelist', choicelist=self.pressureUnits)
        self.setConfigKeyValue(self.defaultSourceName, 'LengthUnits', 'km', valueDataType='choicelist', choicelist=self.lengthUnits)
        self.choiceBoxn2.SetSelection(0)
        self.setCountry(0)
        if self.display_MSLP_tab:
            self.resetMslpMonthList()

    def setCountry(self, event):
        country_selection_number = self.choiceBoxn2.GetSelection()
        if country_selection_number > 0:
            country_selection = self.columnChoices['Countries'][country_selection_number]
            self.choiceBoxn2b.Clear()
            self.choiceBoxn2bb.Clear()
            self.columnChoices['Divisions'] = ['']
            self.columnChoices['Locations'] = ['']
            self.choiceBoxn2bb.Enable(False)
            self.textCtrl3113.SetValue('')
            self.textCtrl3115.SetValue('')
            self.textCtrl3123.SetValue('')
            self.textCtrl3125.SetValue('')
            self.label3111.Enable(False)
            self.label3112.Enable(False)
            self.label3114.Enable(False)
            self.label3121.Enable(False)
            self.label3122.Enable(False)
            self.label3124.Enable(False)
            self.label3111n3.Enable(False)
            self.textCtrl3113.Enable(False)
            self.textCtrl3115.Enable(False)
            self.textCtrl3123.Enable(False)
            self.textCtrl3125.Enable(False)
            self.sqlcur.execute('select placename from localities where parentcountry=? and placetype<>?',(country_selection,'locality'))
            self.columnChoices['Divisions'] = [z[0] for z in self.sqlcur.fetchall()]
            # Remove duplicates
            self.columnChoices['Divisions'] = list(set(self.columnChoices['Divisions']))
            self.columnChoices['Divisions'].sort()
            self.setConfigKeyValue('Region', 'LocalityID', '-99999')
            self.setConfigKeyValue('Region', 'LocalityName', '')
            for divName in self.columnChoices['Divisions']:
                self.choiceBoxn2b.Append(divName)
            self.choiceBoxn2b.Enable(True)
        else:
            self.choiceBoxn2b.Clear()
            self.choiceBoxn2bb.Clear()
            self.columnChoices['Divisions'] = ['']
            self.columnChoices['Locations'] = ['']
            self.choiceBoxn2b.Enable(False)
            self.choiceBoxn2bb.Enable(False)
            self.textCtrl3113.SetValue('')
            self.textCtrl3115.SetValue('')
            self.textCtrl3123.SetValue('')
            self.textCtrl3125.SetValue('')
            self.label3111.Enable(False)
            self.label3112.Enable(False)
            self.label3114.Enable(False)
            self.label3121.Enable(False)
            self.label3122.Enable(False)
            self.label3124.Enable(False)
            self.label3111n3.Enable(False)
            self.textCtrl3113.Enable(False)
            self.textCtrl3115.Enable(False)
            self.textCtrl3123.Enable(False)
            self.textCtrl3125.Enable(False)
            self.setConfigKeyValue('Region', 'LocalityID', '-99999')
            self.setConfigKeyValue('Region', 'LocalityName', '')


    def setDivision(self, event):
        country_selection_number = self.choiceBoxn2.GetSelection()
        division_selection_number = self.choiceBoxn2b.GetSelection()
        if division_selection_number >= 0:
            country_selection = self.columnChoices['Countries'][country_selection_number]
            division_selection = self.columnChoices['Divisions'][division_selection_number]
            self.choiceBoxn2bb.Clear()
            self.columnChoices['Locations'] = ['']
            self.textCtrl3113.SetValue('')
            self.textCtrl3115.SetValue('')
            self.textCtrl3123.SetValue('')
            self.textCtrl3125.SetValue('')
            self.label3111.Enable(False)
            self.label3112.Enable(False)
            self.label3114.Enable(False)
            self.label3121.Enable(False)
            self.label3122.Enable(False)
            self.label3124.Enable(False)
            self.label3111n3.Enable(False)
            self.textCtrl3113.Enable(False)
            self.textCtrl3115.Enable(False)
            self.textCtrl3123.Enable(False)
            self.textCtrl3125.Enable(False)
            self.sqlcur.execute('select placename from localities where parentcountry=? and parentdivision=? and placetype=? order by population desc limit 50',(country_selection,division_selection,'locality'))
            self.columnChoices['Locations'] = [z[0] for z in self.sqlcur.fetchall()]
            self.columnChoices['Locations'].sort()
            self.setConfigKeyValue('Region', 'LocalityID', '-99999')
            self.setConfigKeyValue('Region', 'LocalityName', '')
            for locName in self.columnChoices['Locations']:
                self.choiceBoxn2bb.Append(locName)
            self.choiceBoxn2bb.Enable(True)
        else:
            self.choiceBoxn2bb.Clear()
            self.columnChoices['Locations'] = ['']
            self.choiceBoxn2bb.Enable(False)
            self.textCtrl3113.SetValue('')
            self.textCtrl3115.SetValue('')
            self.textCtrl3123.SetValue('')
            self.textCtrl3125.SetValue('')
            self.label3111.Enable(False)
            self.label3112.Enable(False)
            self.label3114.Enable(False)
            self.label3121.Enable(False)
            self.label3122.Enable(False)
            self.label3124.Enable(False)
            self.label3111n3.Enable(False)
            self.textCtrl3113.Enable(False)
            self.textCtrl3115.Enable(False)
            self.textCtrl3123.Enable(False)
            self.textCtrl3125.Enable(False)
            self.setConfigKeyValue('Region', 'LocalityID', '-99999')
            self.setConfigKeyValue('Region', 'LocalityName', '')
          


    def setLocation(self, event):
        country_selection_number = self.choiceBoxn2.GetSelection()
        division_selection_number = self.choiceBoxn2b.GetSelection()
        location_selection_number = self.choiceBoxn2bb.GetSelection()
        country_selection = self.columnChoices['Countries'][country_selection_number]
        division_selection = self.columnChoices['Divisions'][division_selection_number]
        location_selection = self.columnChoices['Locations'][location_selection_number]
        if location_selection_number >= 0:
            #self.columnChoices['Locations']
            self.sqlcur.execute('select lon, lat, placeID from localities where parentcountry=? and parentdivision=? and placename=? and placetype=?',(country_selection,division_selection,location_selection,'locality'))
            lon, lat, placeID = self.sqlcur.fetchone()
            self.textCtrl3113.SetValue(str(math.floor(lon - 5)))
            self.textCtrl3115.SetValue(str(math.ceil(lon + 5)))
            self.textCtrl3123.SetValue(str(math.floor(lat - 5)))
            self.textCtrl3125.SetValue(str(math.ceil(lat + 5)))
            self.setConfigKeyValue('Region', 'LocalityID', placeID)
            self.setConfigKeyValue('Region', 'LocalityName', location_selection + ', ' + division_selection + ', ' + country_selection + '.')
            self.label3111.Enable(True)
            self.label3112.Enable(True)
            self.label3114.Enable(True)
            self.label3121.Enable(True)
            self.label3122.Enable(True)
            self.label3124.Enable(True)
            self.label3111n3.Enable(True)
            self.textCtrl3113.Enable(True)
            self.textCtrl3115.Enable(True)
            self.textCtrl3123.Enable(True)
            self.textCtrl3125.Enable(True)


    def setPlottingUnits(self, event):
        selection_number = self.choiceBoxn822.GetSelection()
        unit_selection = self.speedUnits[selection_number]
        self.setConfigKeyValue('HazardInterface', 'PlotSpeedUnits', unit_selection)

    def loadSampleData(self):
        # Reset grid values
        for k1 in range(self.gridMaxRows):
            for k2 in range(self.gridMaxCols):
                self.grid_inputCSV.SetCellValue(k1, k2, '')

        # If not a file name, reset columns size
        if len(os.path.basename(self.datafilename)) == 0:
            self.grid_inputCSV.SetSize((self.gridColSize*self.gridMaxCols, self.gridRowSize*self.gridMaxRows))
            self.scrolledWindow122.SetVirtualSize((self.gridColSize*self.gridMaxCols, self.gridRowSize*self.gridMaxRows))
            self.scrolledWindow122.Scroll(0, 0)
        else:
            try:
                datafile = open(self.datafilename, 'rb')
            except IOError:
                self.grid_inputCSV.SetSize((self.gridColSize*self.gridMaxCols, self.gridRowSize*self.gridMaxRows))
                self.scrolledWindow122.SetVirtualSize((self.gridColSize*self.gridMaxCols, self.gridRowSize*self.gridMaxRows))
                self.scrolledWindow122.Scroll(0, 0)
                return

            csvReader = csv.reader(datafile, delimiter=',')

            inputsource = self.getConfigKeyValue('DataProcess', 'Source')
            numheadinglines = int(self.getConfigKeyValue(inputsource, 'NumberOfHeadingLines'))
            for k in range(numheadinglines):
                fileline = csvReader.next()

            max_cols = 0
            for k1 in range(self.gridMaxRows):
                try:
                    fileline = csvReader.next()
                except StopIteration:
                    break
                max_cols = max(len(fileline), max_cols)
                for k2 in range(len(fileline)):
                    self.grid_inputCSV.SetCellValue(k1, k2, fileline[k2].strip())
            self.grid_inputCSV.SetSize((self.gridColSize*self.gridMaxCols, self.gridRowSize*self.gridMaxRows))
            self.scrolledWindow122.SetVirtualSize((self.gridColSize*max_cols, self.gridRowSize*self.gridMaxRows))
            datafile.close()

            # Clear settings for columns that exceed the number of input columns
            # This prevents a bug caused when changing the input file to another with less columns
            colList = self.getConfigKeyValue(inputsource, 'Columns')
            colList_trunc = ','.join(colList.split(',')[:max_cols])
            updatedCols = (',' + colList_trunc).replace(',,',',skip,').replace(',,',',skip,').lstrip(',')
            self.setConfigKeyValue(inputsource, 'Columns', updatedCols)
            self.refreshColumns()
            self.scrolledWindow122.Scroll(0, 0)

    def checkValue(self, event, objId=None):
        if objId is None:
            objId = event.Id

        # Load field settings
        field = self.pointers[self.idMap[objId]][0]
        fieldValue = field.GetValue()
        section = self.pointers[self.idMap[objId]][1]
        key = self.pointers[self.idMap[objId]][2]
        subDictKey = self.pointers[self.idMap[objId]][3]

        # Load corresponding configuration dictionary setting
        keySetting = None
        if not (section == 0):
            if not (key == 0):
                keySetting = self.configDict[section][key]
                settingValue = keySetting[0]
                settingDType = keySetting[1]
                minVal = keySetting[2]
                maxVal = keySetting[3]

        if keySetting is not None:

            if settingDType == 'float' or settingDType == 'int':
                validKeyValue = self.setConfigKeyValue(section, key, fieldValue)

                if not validKeyValue:
                    if settingDType == 'int':
                        numLabel = 'an integer'
                    else:
                        numLabel = 'a number'
                    msg_header = "Value must be " + numLabel
                    field.SetValue(str(settingValue))
                    if (minVal is not None) and (maxVal is not None):
                        self.raiseValueWarning(msg_header + " between " + str(minVal) + " and " + str(maxVal))
                    elif (minVal is not None):
                        self.raiseValueWarning(msg_header + " greater than " + str(minVal))
                    elif (minVal is not None) and (maxVal is not None):
                        self.raiseValueWarning(msg_header + " less than " + str(maxVal))
                    else:
                        self.raiseValueWarning(msg_header, style=wx.OK|wx.ICON_INFORMATION)
            elif settingDType == 'xystep':
                if subDictKey == 'x':
                    newSetting = "{'x':" + fieldValue + ", 'y':" + str(settingValue['y']) + "}"
                elif subDictKey == 'y':
                    newSetting = "{'x':" + str(settingValue['x']) + ", 'y':" + fieldValue + "}"
                validKeyValue = self.setConfigKeyValue(section, key, newSetting)
                if not validKeyValue:
                    field.SetValue(str(settingValue[subDictKey]))
                    self.raiseValueWarning('Value must be a number greater than 0')
            elif settingDType == 'domain':
                xMinStr = "{'xMin':" + str(settingValue['xMin']) + ","
                xMaxStr = "'xMax':" + str(settingValue['xMax']) + ","
                yMinStr = "'yMin':" + str(settingValue['yMin']) + ","
                yMaxStr = "'yMax':" + str(settingValue['yMax']) + "}"

                if subDictKey == 'xMin':
                    newSetting = "{'xMin':" + fieldValue + "," + xMaxStr + yMinStr + yMaxStr
                elif subDictKey == 'xMax':
                    newSetting = xMinStr + "'xMax':" + fieldValue + "," + yMinStr + yMaxStr
                elif subDictKey == 'yMin':
                    newSetting = xMinStr + xMaxStr + "'yMin':" + fieldValue + "," + yMaxStr
                elif subDictKey == 'yMax':
                    newSetting = xMinStr + xMaxStr + yMinStr + "'yMax':" + fieldValue + "}"

                validKeyValue = self.setConfigKeyValue(section, key, newSetting)
                if not validKeyValue:
                    field.SetValue(str(settingValue[subDictKey]))
                    if subDictKey == 'xMin' or subDictKey == 'xMax':
                        self.raiseValueWarning("Value must be a number between 0 and 360")
                    elif subDictKey == 'yMin' or subDictKey == 'yMax':
                        self.raiseValueWarning("Value must be a number between -90 and 90")
            elif settingDType == 'intlist':
                validKeyValue = self.setConfigKeyValue(section, key, fieldValue)
                settingValue = self.getConfigKeyValue(section, key)
                field.SetValue(settingValue.strip('[]'))
                if not validKeyValue:
                    self.raiseValueWarning("A list of integers between " + str(minVal) + " and " + str(maxVal) + " is required")

    def raiseValueWarning(self, string, style=wx.OK|wx.ICON_INFORMATION):
        self.valueWarningRaised = True
        MsgDlg(self, string, "Warning", style)

    def refreshGUIsettings( self ):
        # Input Settings
        inputsource = self.getConfigKeyValue('DataProcess', 'Source')

        self.refreshColumns()
        numHeadingLines = self.getConfigKeyValue(inputsource, 'NumberOfHeadingLines')
        self.m_spinCtrl1.SetValue(int(numHeadingLines))

        self.LoadCSVFile()

        # Output Settings
        self.textCtrl212.SetValue(self.getConfigKeyValue('Output', 'Path'))

        # Domain Settings
        self.textCtrl3113.SetValue(self.getConfigKeyValue('Region', 'gridLimit', 'xMin'))
        self.textCtrl3115.SetValue(self.getConfigKeyValue('Region', 'gridLimit', 'xMax'))
        self.textCtrl3123.SetValue(self.getConfigKeyValue('Region', 'gridLimit', 'yMin'))
        self.textCtrl3125.SetValue(self.getConfigKeyValue('Region', 'gridLimit', 'yMax'))

        # RMW Settings
        if self.display_RMW_tab:
            self.textCtrl9123.SetValue(self.getConfigKeyValue('RMW', 'mean'))
            self.textCtrl9125.SetValue(self.getConfigKeyValue('RMW', 'sigma'))

            if self.getConfigKeyValue('RMW', 'GetRMWDistFromInputData'):
                self.radioBtn913.SetValue(True)
                self.onRMWsel2(None)
            else:
                self.radioBtn911.SetValue(True)
                self.onRMWsel1(None)

        # MSLP Settings
        if self.display_MSLP_tab:
            if self.yearlyAverageChecked:
                self.radioBtn411.SetValue(True)
                self.label41321.Enable(False)
                for k in range(12):
                    self.month_sel[k].Enable(False)
            else:
                for k in range(12):
                    self.month_sel[k].Enable(True)
                self.radioBtn412.SetValue(True)
                self.label41321.Enable(True)

        # Stats Interface Settings
        if self.display_stats_tab:
            kdeType1D = self.getConfigKeyValue('StatInterface', 'kdeType')
            if kdeType1D == self.KDETypes[0]:
                self.radioBtn5112.SetValue(True)
            elif kdeType1D == self.KDETypes[1]:
                self.radioBtn5113.SetValue(True)
            elif kdeType1D == self.KDETypes[2]:
                self.radioBtn5114.SetValue(True)
            elif kdeType1D == self.KDETypes[3]:
                self.radioBtn5115.SetValue(True)

            kdeType2D = self.getConfigKeyValue('StatInterface', 'kde2DType')
            if kdeType2D == self.KDETypes[0]:
                self.radioBtn5122.SetValue(True)
            elif kdeType2D == self.KDETypes[1]:
                self.radioBtn5123.SetValue(True)
            elif kdeType2D == self.KDETypes[2]:
                self.radioBtn5124.SetValue(True)
            elif kdeType2D == self.KDETypes[3]:
                self.radioBtn5125.SetValue(True)

            self.textCtrl5213.SetValue(self.getConfigKeyValue('StatInterface', 'gridSpace', 'x'))
            self.textCtrl5215.SetValue(self.getConfigKeyValue('StatInterface', 'gridSpace', 'y'))
            self.textCtrl5223.SetValue(self.getConfigKeyValue('StatInterface', 'gridInc', 'x'))
            self.textCtrl5225.SetValue(self.getConfigKeyValue('StatInterface', 'gridInc', 'y'))

        # Track Generator Settings
        self.textCtrl6112.SetValue(self.getConfigKeyValue('TrackGenerator', 'NumSimulations'))
        #self.textCtrl6122.SetValue(self.getConfigKeyValue('TrackGenerator', 'YearsPerSimulation'))
        #self.textCtrl6212.SetValue(self.getConfigKeyValue('TrackGenerator', 'Frequency'))

        # Windfield Settings
        if self.display_windfield_tab:
            profile_setting = self.getConfigKeyValue('WindfieldInterface', 'profileType')
            if profile_setting == self.profileTypes[0]:
                self.radioBtn7111.SetValue(True)
                self.onOtherProfile(None)
            elif profile_setting == self.profileTypes[1]:
                self.radioBtn7112.SetValue(True)
                self.onDoubleHollandProfile(None)
            elif profile_setting == self.profileTypes[2]:
                self.radioBtn7113.SetValue(True)
                self.onHollandProfile(None)
            elif profile_setting == self.profileTypes[3]:
                self.radioBtn7114.SetValue(True)
                self.onOtherProfile(None)
            elif profile_setting == self.profileTypes[4]:
                self.radioBtn7115.SetValue(True)
                self.onOtherProfile(None)
            elif profile_setting == self.profileTypes[5]:
                self.radioBtn7116.SetValue(True)
                self.onOtherProfile(None)

            self.textCtrl7122.SetValue(self.getConfigKeyValue('WindfieldInterface', 'beta'))
            #self.textCtrl7124.SetValue(self.getConfigKeyValue('WindfieldInterface', 'beta1'))
            #self.textCtrl7126.SetValue(self.getConfigKeyValue('WindfieldInterface', 'beta2'))

            #boundarylayer_setting = self.getConfigKeyValue('WindfieldInterface', 'windFieldType')
            #if boundarylayer_setting == self.windFieldTypes[0]:
            #    self.radioBtn721.SetValue(True)
            #    self.onKepert(None)
            #elif boundarylayer_setting == self.windFieldTypes[1]:
            #    self.radioBtn722.SetValue(True)
            #    self.onMcConochieOrHubbert(None)
            #elif boundarylayer_setting == self.windFieldTypes[2]:
            #    self.radioBtn723.SetValue(True)
            #    self.onMcConochieOrHubbert(None)

            #self.textCtrl725.SetValue(self.getConfigKeyValue('WindfieldInterface', 'thetaMax'))
            #self.textCtrl732.SetValue(self.getConfigKeyValue('WindfieldInterface', 'Resolution'))

        # Hazard Settings
        self.textCtrl812.SetValue(self.getConfigKeyValue('HazardInterface','Years').strip('[]'))
        plotUnitSelection = self.getConfigKeyValue('HazardInterface', 'PlotSpeedUnits')
        plotUnitSelectionIndex = self.speedUnits.index(plotUnitSelection)
        self.choiceBoxn822.SetSelection(int(plotUnitSelectionIndex))

    def refreshColumns(self):
        inputsource = self.getConfigKeyValue('DataProcess', 'Source')
        columnsStr = self.getConfigKeyValue(inputsource, 'Columns')
        columnsList = [k.strip() for k in columnsStr.split(',')]
        for k in range(len(self.colmap)):
            if k < len(columnsList):
                selectionNo = numpy.flatnonzero(numpy.array([n == columnsList[k] for n in self.columnChoices['refName']]))
                if len(selectionNo) == 0:
                    selectionNo = 0
                else:
                    selectionNo = selectionNo[0]
            else:
                selectionNo = 0
            self.colmap[k].SetSelection(selectionNo)

    def menu_sel_file_new( self, event ):
        warningRaised = self.checkAllValues()
        if not warningRaised:
            self.createNewConfigFile(None)

    def menu_sel_file_open( self, event ):
        warningRaised = self.checkAllValues()
        if not warningRaised:
            self.LoadConfigFile(None)

    def menu_sel_file_save( self, event, newfilename=None ):
        warningRaised = self.checkAllValues()
        if not warningRaised:

            if newfilename is None:
                configFilename = self.textCtrl012.GetValue()
                if configFilename == self.defaultconfigname:
                    return self.menu_sel_file_saveas(0)
            else:
                configFilename = newfilename

            self.updateColumnsKeyValue()

            # Copy settings across modules
            self.setConfigKeyValue('TrackGenerator', 'gridSpace', self.getConfigKeyValue('StatInterface', 'gridSpace'))
            self.setConfigKeyValue('TrackGenerator', 'gridInc', self.getConfigKeyValue('StatInterface', 'gridInc'))
            self.setConfigKeyValue('HazardInterface', 'Resolution', self.getConfigKeyValue('WindfieldInterface', 'Resolution'))
            self.setConfigKeyValue('WindfieldInterface', 'NumberofFiles', self.getConfigKeyValue('TrackGenerator', 'NumSimulations'))
            self.setConfigKeyValue('HazardInterface', 'NumSim', self.getConfigKeyValue('TrackGenerator', 'NumSimulations'))
            self.setConfigKeyValue('HazardInterface', 'YearsPerSimulation', self.getConfigKeyValue('TrackGenerator', 'YearsPerSimulation'))
            #self.setConfigKeyValue('StatInterface', 'gridLimit', self.getConfigKeyValue('TrackGenerator', 'gridLimit'))

            if self.display_MSLP_tab:
                if self.yearlyAverageChecked:
                    self.setConfigKeyValue('Input', 'MSLPGrid', '[1 2 3 4 5 6 7 8 9 10 11 12]')
                else:
                    monthlist = [(k+1) for k in range(12) if self.month_sel[k].GetValue()]
                    validKeyValue = self.setConfigKeyValue('Input', 'MSLPGrid', str(monthlist))
                    if not validKeyValue:
                        self.setConfigKeyValue('Input', 'MSLPGrid', '[1 2 3 4 5 6 7 8 9 10 11 12]')

            try:
                f = open(configFilename, 'w')
                for section in self.configDictSections:
                    f.write('[' + section + ']' + '\n')
                    settings = self.configDict[section].keys()
                    recOrder = [self.configDict[section][setting][5] for setting in settings]
                    recOrderSorted = [numpy.flatnonzero(k==recOrder)[0] for k in numpy.sort(recOrder)]
                    settings = [settings[k] for k in recOrderSorted]
                    for setting in settings:
                        settingValue = deepcopy(self.configDict[section][setting][0])
                        if self.configDict[section][setting][1] == 'domain':
                            settingValue = self.convertDictKeyValues2Strings(settingValue)
                            outputline = setting + '=' + '{' + "'xMin':" + settingValue['xMin'] + \
                                        ",'xMax':" + settingValue['xMax'] + ",'yMin':" + \
                                        settingValue['yMin'] + ",'yMax':" + settingValue['yMax'] + '}'
                            f.write(outputline + '\n')
                        elif self.configDict[section][setting][1] == 'xystep':
                            settingValue = self.convertDictKeyValues2Strings(settingValue)
                            outputline = setting + '=' + '{' + "'x':" + settingValue['x'] + ",'y':" + \
                                        settingValue['y'] + '}'
                            f.write(outputline + '\n')
                        elif self.configDict[section][setting][1] == 'intlist':
                            outputline = setting + '=' + ','.join([str(k) for k in settingValue])
                            f.write(outputline + '\n')
                        else:
                            # Handle unicode place names by converting to ASCII
                            if isinstance(settingValue, unicode):
                                settingValue = unicodedata.normalize('NFKD', settingValue).encode('ascii', 'ignore')
                            
                            settingValue = str(settingValue)                            
                            if setting[0] == ';':   # Check if comment line or setting
                                outputline = setting
                            else:
                                if not (settingValue.find('{output directory}')==-1):
                                    outputPath = self.getConfigKeyValue('Output', 'Path')
                                    settingValue = settingValue.replace('{output directory}', outputPath)
                                    settingValue = os.path.abspath(settingValue)
                                if not (settingValue.find('{tcrm directory}')==-1):
                                    settingValue = settingValue.replace('{tcrm directory}', self.tcrmPath)
                                    settingValue = os.path.abspath(settingValue)
                                if not (settingValue.find('{filename}')==-1):
                                    basefilename = os.path.basename(configFilename).rstrip('.ini')
                                    settingValue = settingValue.replace('{filename}', basefilename)
                                    settingValue = os.path.abspath(settingValue)
                                outputline = setting + '=' + settingValue
                            f.write(outputline + '\n')
                    f.write('\n')
                f.close()
                self.settingsChanged = False
                self.refreshGUIsettings()
                self.SetStatusText(" File saved")
                time.sleep(0.4)
                self.SetStatusText("")
            except KeyError:
                MsgDlg(self, "An error occured when writing to the file!   ", caption='Warning', style=wx.OK|wx.ICON_ERROR)
            except IOError, e:
                if e.errno == 13:
                    MsgDlg(self, "You do not have permission to save in this directory.", caption='Warning', style=wx.OK|wx.ICON_ERROR)
                    return True
                else:
                    MsgDlg(self, "An error occured when writing to the file!   ", caption='Warning', style=wx.OK|wx.ICON_ERROR)

    def menu_sel_file_saveas( self, event ):
        warningRaised = self.checkAllValues()
        if not warningRaised:

            filename = self.textCtrl012.GetValue()
            if filename == self.defaultconfigname:
                filename = '*.ini'
            dlg = wx.FileDialog(self, 'Save As', '.', filename, 'Configuration Files (*.ini)|*.ini', wx.SAVE | wx.OVERWRITE_PROMPT)

            vetoExit = False

            while 1:
                if dlg.ShowModal() == wx.ID_OK:
                    newfilename = dlg.GetPath()
                    noWritePermission = self.menu_sel_file_save( 0, newfilename)
                    if noWritePermission:
                        # If no write permission, then re-open SaveAs dialog
                        pass
                    else:
                        self.textCtrl012.SetValue(newfilename)
                        self.backpanel.Update()                    
                        break
                else:
                    # Flag to prevent program exit if user cancelled file save
                    vetoExit = True
                    break

            dlg.Destroy()
            return vetoExit

    def menu_sel_file_exit( self, event ):
        self.Close()

    def save_and_continue( self, event ):
        self.menu_sel_file_save(0)
        tab_number = self.notebook_panel.GetSelection()
        if tab_number < self.notebook_panel.GetPageCount() - 1:
            self.notebook_panel.SetSelection(tab_number + 1)

    def createNewConfigFile( self, event ):
        if self.settingsChanged:
            filename = self.textCtrl012.GetValue()
            if len(filename) > 0:
                selFlag = MsgDlg(self, 'Do you want to save the changes to ' + filename + '?', caption='TCRM Configuration Editor', style=wx.YES_NO|wx.CANCEL|wx.ICON_EXCLAMATION)
                if selFlag == wx.ID_YES:
                    exitFlag = self.menu_sel_file_save(0)
                    if exitFlag is True:
                        return
                elif selFlag == wx.ID_NO:
                    pass
                else:
                    return
        self.textCtrl012.SetValue(self.defaultconfigname)
        self.initialiseSettings()
        self.refreshGUIsettings()
        self.notebook_panel.Enable( True )
        self.tabsEnabled = True
        self.file_save.Enable( True )
        self.file_saveas.Enable( True )
        self.settingsChanged = False
        self.m_button8.Enable(True)
        self.m_button9.Enable(True)
        self.notebook_panel.SetSelection(0)
        self.updateUserInfo()

    def onMouseOverButton013( self, event ):
        self.SetStatusText(u" Create a new configuration file")

    def onMouseOverButton014( self, event ):
        self.SetStatusText(u" Open an existing configuration file")

    def onMouseLeave( self, event ):
        self.SetStatusText("")

    def headerLinesSettingChange( self, event ):
        inputsource = self.getConfigKeyValue('DataProcess', 'Source')
        self.setConfigKeyValue(inputsource, 'NumberOfHeadingLines', str(self.m_spinCtrl1.GetValue()))
        self.onSettingsChanged()
        # Refresh sample data
        self.loadSampleData()

    def updateColumnsKeyValue(self):
        inputsource = self.getConfigKeyValue('DataProcess', 'Source')            
        # Truncates unassigned columns on right and converts rest into skip value
        colList = self.getConfigKeyValue(inputsource, 'Columns')
        updatedCols = (',' + colList).replace(',,',',skip,').replace(',,',',skip,').lstrip(',')
        self.setConfigKeyValue(inputsource, 'Columns', updatedCols)
        
    def LoadConfigFile( self, event ):
        if self.settingsChanged:
            filename = self.textCtrl012.GetValue()
            if len(filename) > 0:
                selFlag = MsgDlg(self, 'Do you want to save the changes to ' + filename + '?', caption='TCRM Configuration Editor', style=wx.YES_NO|wx.CANCEL|wx.ICON_EXCLAMATION)
                if selFlag == wx.ID_YES:
                    exitFlag = self.menu_sel_file_save(0)
                    if exitFlag is True:
                        return
                elif selFlag == wx.ID_NO:
                    pass
                else:
                    return
        dlg = wx.FileDialog(self, 'Open', '.', '', 'Configuration Files (*.ini)|*.ini', wx.OPEN | wx.FILE_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            fullfilename = dlg.GetPath()
            self.textCtrl012.SetValue(fullfilename)
            self.readConfigFile(fullfilename)
            self.notebook_panel.Enable( True )
            self.tabsEnabled = True
            self.file_save.Enable( True )
            self.file_saveas.Enable( True )
            self.settingsChanged = False
            self.m_button8.Enable(True)
            self.m_button9.Enable(True)
            self.notebook_panel.SetSelection(0)
            self.updateUserInfo()
        dlg.Destroy()

    def LoadDirPicker( self, event ):
        dlg = wx.DirDialog(self, 'Select TCRM Output Folder', self.getConfigKeyValue('Output', 'Path'))
        if dlg.ShowModal() == wx.ID_OK:
            self.setConfigKeyValue('Output', 'Path', dlg.GetPath())
            self.textCtrl212.SetValue(self.getConfigKeyValue('Output', 'Path'))
        dlg.Destroy()

    def LoadCSVFile(self):
        self.datafilename = self.getConfigKeyValue('DataProcess', 'InputFile')
        # If input file has no path information, default to tcrm input folder
        if len(os.path.dirname(self.datafilename)) == 0:
            self.datafilename = os.path.join(self.tcrmInputDir, self.datafilename)
        # Display text if filename, otherwise suppress text (since it will represent only the default directory)
        if len(os.path.basename(self.datafilename)) > 0:
            self.textCtrl1119.SetValue(self.datafilename)
            self.loadSampleData()
        else:
            self.textCtrl1119.SetValue('')
            self.loadSampleData()

    def OnProjectExit(self, event):
        warningRaised = self.checkAllValues()
        if warningRaised:
            event.Veto();
            return
        if self.settingsChanged:
            filename = self.textCtrl012.GetValue()
            if event.CanVeto() & (len(filename) > 0):
                selFlag = MsgDlg(self, 'Do you want to save the changes to ' + filename + '?', caption='TCRM Configuration Editor', style=wx.YES_NO|wx.CANCEL|wx.ICON_EXCLAMATION)
                if selFlag == wx.ID_YES:
                    exitFlag = self.menu_sel_file_save(0)
                    if exitFlag is None:
                        pass
                    elif exitFlag is True:
                        # User cancelled file save => veto program exit
                        event.Veto();
                        return
                elif selFlag == wx.ID_NO:
                    pass
                else:
                    event.Veto();
                    return
        self.sqlcur.close()
        self.sqlcon.close()
        event.Skip()

    def SelectCSVFile(self, event):
        current_filename = self.getConfigKeyValue('DataProcess', 'InputFile')
        dirname, filename = os.path.split(current_filename)
        selected_filename = wx.FileSelector(u"Select an input CSV (comma delimited) track file", dirname, filename, wx.EmptyString, "*.csv", wx.FD_FILE_MUST_EXIST) #|wx.FLP_OPEN)
        if len(selected_filename):
            self.setConfigKeyValue('DataProcess', 'InputFile', selected_filename)
            new_filename = self.getConfigKeyValue('DataProcess', 'InputFile')
            # Test if file has changed
            if new_filename != current_filename:
                self.LoadCSVFile()                

    def OnSeasonalAverage(self, event):
        self.onSettingsChanged()
        self.yearlyAverageChecked = False
        self.refreshGUIsettings()

    def OnYearlyAverage(self, event):
        self.onSettingsChanged()
        self.yearlyAverageChecked = True
        self.refreshGUIsettings()

    def onHollandProfile(self, event):
        self.updateProfileSetting()
        self.label7121.Enable( True )
        self.textCtrl7122.Enable( True )
        #self.label7123.Enable( False )
        #self.textCtrl7124.Enable( False )
        #self.label7125.Enable( False )
        #self.textCtrl7126.Enable( False )

    def onDoubleHollandProfile(self, event):
        self.updateProfileSetting()
        self.label7121.Enable( False )
        self.textCtrl7122.Enable( False )
        #self.label7123.Enable( True )
        #self.textCtrl7124.Enable( True )
        #self.label7125.Enable( True )
        #self.textCtrl7126.Enable( True )

    def onOtherProfile(self, event):
        self.updateProfileSetting()
        self.label7121.Enable( False )
        self.textCtrl7122.Enable( False )
        #self.label7123.Enable( False )
        #self.textCtrl7124.Enable( False )
        #self.label7125.Enable( False )
        #self.textCtrl7126.Enable( False )

    def updateProfileSetting(self):
        if self.radioBtn7111.GetValue():
            profileselection = 0
        elif self.radioBtn7112.GetValue():
            profileselection = 1
        elif self.radioBtn7113.GetValue():
            profileselection = 2
        elif self.radioBtn7114.GetValue():
            profileselection = 3
        elif self.radioBtn7115.GetValue():
            profileselection = 4
        elif self.radioBtn7116.GetValue():
            profileselection = 5
        self.setConfigKeyValue('WindfieldInterface', 'profileType', self.profileTypes[profileselection])

    def updateBoundaryLayerSetting(self):
        if self.radioBtn721.GetValue():
            blselection = 0
        elif self.radioBtn722.GetValue():
            blselection = 1
        elif self.radioBtn723.GetValue():
            blselection = 2
        self.setConfigKeyValue('WindfieldInterface', 'windFieldType', self.windFieldTypes[blselection])

    def updateKDE1DSettings(self, event):
        if self.radioBtn5112.GetValue():
            kde1Dselection = 0
        elif self.radioBtn5113.GetValue():
            kde1Dselection = 1
        elif self.radioBtn5114.GetValue():
            kde1Dselection = 2
        elif self.radioBtn5115.GetValue():
            kde1Dselection = 3
        self.setConfigKeyValue('StatInterface', 'kdeType', self.KDETypes[kde1Dselection])

    def updateKDE2DSettings(self, event):
        if self.radioBtn5122.GetValue():
            kde2Dselection = 0
        elif self.radioBtn5123.GetValue():
            kde2Dselection = 1
        elif self.radioBtn5124.GetValue():
            kde2Dselection = 2
        elif self.radioBtn5125.GetValue():
            kde2Dselection = 3
        self.setConfigKeyValue('StatInterface', 'kde2DType', self.KDETypes[kde2Dselection])

    def onTabChange(self, event):
        if self.tabsEnabled:
            warningRaised = self.checkAllValues()
            if warningRaised:
                event.Veto();
                return

    def onTabChanged(self, event):
        self.updateUserInfo()

    def updateUserInfo(self):
        self.m_staticText36.SetLabel(self.userhelptext[self.notebook_panel.GetSelection()])
        self.bSizer421.FitInside( self.m_scrolledWindow2 )

    def onColMapChange(self, event):
        if not self.columnDuplicates():
            colsettings = [self.colmap[k].GetSelection() for k in range(self.gridMaxCols)]
            colnamesettings = [self.columnChoices['refName'][colsettings[k]] for k in range(self.gridMaxCols)]
            inputsource = self.getConfigKeyValue('DataProcess', 'Source')
            self.setConfigKeyValue(inputsource, 'Columns', ','.join(colnamesettings).rstrip(','), valueDataType='str')

            selIdx = self.colmap[event.Id - 9000].GetSelection()
            changedcolvalue = self.columnChoices['refName'][selIdx]
            unitkey = self.columnChoices['unitMapping'][selIdx]

            if unitkey == '':
                pass
            elif unitkey == 'DateFormat':
                dlg = DateFormatDialog(self)
                result = dlg.ShowModal()
                dlg.Destroy()
            else:
                selectionOptions = self.configDict[inputsource][unitkey][4]
                inputUnits = self.getConfigKeyValue(inputsource, unitkey)
                dlg = UnitSelectionBox(self, unitkey, selectionOptions, inputUnits)
                result = dlg.ShowModal()
                dlg.Destroy()
                if result in range(len(selectionOptions)):
                    self.setConfigKeyValue(inputsource, unitkey, selectionOptions[result])
        else:
            self.raiseValueWarning('Duplicate columns are not allowed', style=wx.OK|wx.ICON_INFORMATION)
            self.refreshColumns()

    def columnDuplicates(self):
        colsettings = deepcopy([self.colmap[k].GetSelection() for k in range(self.gridMaxCols)])
        colnamesettings = [self.columnChoices['refName'][colsettings[k]] for k in range(self.gridMaxCols)]

        columnIdxs = deepcopy(colsettings)
        for k in range(self.gridMaxCols):
            if colnamesettings[k] == 'skip':
                # Allow multiple skip columns
                columnIdxs[k] = 0
        isduplicates = not (sum(columnIdxs) == sum(numpy.unique(columnIdxs)))
        return isduplicates

    def onKepert(self, event):
        self.updateBoundaryLayerSetting()
        #self.label724.Enable( False )
        #self.textCtrl725.Enable( False )

    def onMcConochieOrHubbert(self, event):
        self.updateBoundaryLayerSetting()
        self.label724.Enable( True )
        self.textCtrl725.Enable( True )

    def onRMWsel1(self, event):
        self.setConfigKeyValue('RMW', 'GetRMWDistFromInputData', 'False')
        self.label9122.Enable(True)
        self.textCtrl9123.Enable(True)
        self.label9124.Enable(True)
        self.textCtrl9125.Enable(True)

    def onRMWsel2(self, event):
        self.setConfigKeyValue('RMW', 'GetRMWDistFromInputData', 'True')
        self.label9122.Enable(False)
        self.textCtrl9123.Enable(False)
        self.label9124.Enable(False)
        self.textCtrl9125.Enable(False)

    def onMSLPcheckbox(self, event):
        self.onSettingsChanged()

    def onSettingsChanged(self):
        if not self.settingsChanged:
            self.settingsChanged = True

    def readConfigFile(self, configFilename):
        cp = ConfigParser.ConfigParser()
        cp.optionxform = str

        try:
            configFile = open(configFilename)
            cp.readfp(configFile)
            configFile.close()
        except:
            MsgDlg(self, "Error loading config file", caption='Warning', style=wx.OK|wx.ICON_ERROR)

        self.initialiseSettings()

        for section in cp.sections():
            for key in cp.options(section):
                self.setConfigKeyValue(section, key, cp.get(section, key))

        sourceName = self.getConfigKeyValue('DataProcess', 'Source')

        if not (sourceName == self.defaultSourceName):
            if self.configDict.has_key(sourceName):
                sourceKeys = self.configDict[sourceName].keys()
                for key in sourceKeys:
                    self.setConfigKeyValue(self.defaultSourceName, key, str(self.getConfigKeyValue(sourceName, key)))
                self.configDict.pop(sourceName)
                self.configDict[sourceName] = self.configDict.pop(self.defaultSourceName)
                self.configDictSections.remove(self.defaultSourceName)

        if self.display_MSLP_tab:
            self.resetMslpMonthList()

        # Reset column settings if duplicate columns are listed in config file to avoid introducing bug
        self.refreshGUIsettings()
        if self.columnDuplicates():
            inputsource = self.getConfigKeyValue('DataProcess', 'Source')
            self.setConfigKeyValue(inputsource, 'Columns', '')
            self.refreshColumns()
        placeID = self.getConfigKeyValue('Region', 'LocalityID')
        
        #country_selection_number = self.choiceBoxn2.GetSelection()
        self.sqlcur.execute('select placename, parentcountry, parentdivision from localities where placeID=?', (placeID,))
        localityMatch = self.sqlcur.fetchone()

        if localityMatch is not None:
            placeName, parentCountry, parentDivision = localityMatch
            if parentCountry in self.columnChoices['Countries']:
                self.choiceBoxn2.SetSelection(self.columnChoices['Countries'].index(parentCountry))
                self.setCountry(0)
                self.choiceBoxn2b.SetSelection(self.columnChoices['Divisions'].index(parentDivision))
                self.setDivision(0)
                self.choiceBoxn2bb.SetSelection(self.columnChoices['Locations'].index(placeName))
                self.setLocation(0)
        else:
            self.setConfigKeyValue('Region', 'LocalityID', '-99999')

            
    def resetMslpMonthList(self):
        monthlist = set(self.convertString2IntList(self.getConfigKeyValue('Input', 'MSLPGrid')))
        allmonths = set(range(1,13))
        if monthlist == allmonths:
            self.yearlyAverageChecked = True
        else:
            self.yearlyAverageChecked = False
        for k in monthlist:
            self.month_sel[k - 1].SetValue(True)
        for k in allmonths.difference(monthlist):
            self.month_sel[k - 1].SetValue(False)

    def getConfigKeyValue(self, section, key, subkey=None):
        if subkey is None:
            if self.configDict[section][key][1] == 'logical':
                return self.configDict[section][key][0]
            else:
                return str(self.configDict[section][key][0])
        else:
            return str(self.configDict[section][key][0][subkey])

    def setConfigKeyValue(self, section, key, keyValue, valueDataType=None, minVal=None, maxVal=None, choicelist=[], allowOverwrite=True):
        # Create config section if doesn't exist
        if not self.configDict.has_key(section):
            self.configDict[section] = {}
            self.configDictSections.append(section)

        # If key is pre-initialised, load key settings (unless overwrite settings are given as input)
        keyRecNo = None
        if self.configDict[section].has_key(key):
            currentKeySettings = self.configDict[section][key]
            if valueDataType is None:
                valueDataType = currentKeySettings[1]
            if minVal is None:
                minVal = currentKeySettings[2]
            if maxVal is None:
                maxVal = currentKeySettings[3]
            if choicelist == []:
                choicelist = currentKeySettings[4]
            keyRecNo = currentKeySettings[5]
            allowOverwrite = currentKeySettings[6]

        if keyValue is None:
            keyValue = ''
        elif len(keyValue) == 0:
            pass
        else:
            if valueDataType=='str':
                keyValue = str(keyValue)
            elif (valueDataType=='int') or (valueDataType=='float'):
                try:
                    keyValue = float(keyValue)
                    if valueDataType=='int':
                        keyValue = int(keyValue)
                    if minVal is not None:
                        if keyValue < minVal:
                            keyValue = None
                    if maxVal is not None:
                        if keyValue > maxVal:
                            keyValue = None
                except ValueError:
                    keyValue = None
            elif valueDataType=='logical':
                if keyValue == 'True':
                    keyValue = True
                elif keyValue == 'False':
                    keyValue = False
                else:
                    keyValue = None
            elif valueDataType=='intlist':
                keyValue = self.convertString2IntList(keyValue)
                if keyValue is not None:
                    if any([k < minVal for k in keyValue]) and (minVal is not None):
                        keyValue = None
                    if any([k > maxVal for k in keyValue]) and (maxVal is not None):
                        keyValue = None

        if keyValue is not None:
            if valueDataType=='xystep':
                keyValue = self.convertString2Dict(keyValue)
                if (keyValue.has_key('x') and keyValue.has_key('y')):
                    if ((keyValue['x'] > 0) or (keyValue['x'] == '')) and \
                       ((keyValue['y'] > 0) or (keyValue['y'] == '')):
                        pass
                    else:
                        keyValue = None
                else:
                    keyValue = None
            elif valueDataType=='domain':
                keyValue = self.convertString2Dict(keyValue)
                if (keyValue.has_key('xMin') and keyValue.has_key('xMax') and
                    keyValue.has_key('yMin') and keyValue.has_key('yMax')):
                    if ((keyValue['xMin'] >= 0 and keyValue['xMin'] <= 360) or (keyValue['xMin'] == '')) and \
                       ((keyValue['xMax'] >= 0 and keyValue['xMax'] <= 360) or (keyValue['xMax'] == '')) and \
                       ((keyValue['yMin'] >= -90 and keyValue['yMin'] <= 90) or (keyValue['yMin'] == '')) and \
                       ((keyValue['yMax'] >= -90 and keyValue['yMax'] <= 90) or (keyValue['yMax'] == '')):
                        pass
                    else:
                        keyValue = None
                else:
                    keyValue = None
            elif valueDataType=='choicelist':
                keyValue = str(keyValue)
                if not any([k == keyValue for k in choicelist]):
                    keyValue = None

        # If valid keyValue => add/update TCRM configuration dictionary
        validKeyValue = False
        if keyValue is not None:
            if keyRecNo is None:
                keyRecNo = self.nextRecNo
                self.nextRecNo = self.nextRecNo + 1
                newRecord = True
            else:
                newRecord = False
            if self.configDict[section].has_key(key):
                if not (deepcopy(self.configDict[section][key][0]) == deepcopy(keyValue)):
                    self.onSettingsChanged()
                else:
                    validKeyValue = True
            if newRecord or allowOverwrite:
                self.configDict[section][key] = deepcopy([keyValue, valueDataType, minVal, maxVal, choicelist, keyRecNo, allowOverwrite])
                validKeyValue = True
        return validKeyValue

    def checkAllValues(self):
        self.valueWarningRaised = False
        for k in self.idMap.keys():
            self.checkValue(None, objId=k)
        return self.valueWarningRaised

    def loadHelp(self, event):
        if event.Id == 1010:
            helpText = "Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)\n" + \
                       "Copyright (C) 2011  Geoscience Australia\n\n" + \
                       "This program is free software: you can redistribute it and/or modify " + \
                       "it under the terms of the GNU General Public License as published by " + \
                       "the Free Software Foundation, either version 3 of the License, or " + \
                       "(at your option) any later version." + \
                       "\n\n" + \
                       "This program is distributed in the hope that it will be useful, " + \
                       "but WITHOUT ANY WARRANTY; without even the implied warranty of " + \
                       "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the " + \
                       "GNU General Public License for more details." + \
                       "\n\n" + \
                       "You should have received a copy of the GNU General Public License " + \
                       "along with this program.  If not, see <http://www.gnu.org/licenses/>."
                       
        elif event.Id == 3232:
            helpText = "TCRM requires the user to select two domains: \n\n" + \
                       "The 'windfield/hazard' domain specifies the region to calculate the return-period wind speeds.\n\n" + \
                       "The 'track generator' domain specifies the larger region that encompasses all tracks that enter the windfield/hazard domain. " + \
                       "This domain can be auto-calculated from your input track file by pressing the 'autocalculate' button."
        elif event.Id == 914:
            helpText = "TCRM provides two ways to define the Radius of Maximum Wind (RMW) distribution: \n\n" + \
                       "Either:\n\n" + \
                       "(1) A log-normal distribution can be selected with mean (in kilometres) and sigma values chosen by the user\n\n" + \
                       "or\n\n" + \
                       "(2) The distribution can be calculated from the input data using kernal density estimation.  Note: this requires a RMW column in the input file."
        elif event.Id == 414:
            helpText = "A Mean Sea Level Pressure (MSLP) field is required for producing parametric wind profiles.\n\n" + \
                       "Select either:\n\n" + \
                       "(1) A yearly average climatology\n" + \
                       "or\n" + \
                       "(2) A seasonal climatology.  Choose months corresponding to the tropical cyclone active period.\n\n" + \
                       "Climatology derived from NCEP Reanalysis-2 with date range: 1980-2007."
        elif event.Id == 6222:
            helpText = "The annual frequency of storms must be specified for the 'track generator' domain.\n\n" + \
                       "This can be selected manually or calculated automatically from the input data by pressing the 'autocalculate' button."
        else:
            helpText = "{Missing Text}"

        if event.Id == 1010:
            dlg = HelpPopUp(self, helpText, "About")
        else:
            dlg = HelpPopUp(self, helpText, "Help")
        result = dlg.ShowModal()
        dlg.Destroy()

    def convertString2Dict(self, strInp):
        # Converts a string to a dictionary (avoiding unsafe methods such as eval())
        # Assumes simple {string:float, ...} dictionary format
        try:
            keyValuePairs = [k.split(':') for k in strInp.strip('{}').split(',')]
            keys = [str(k[0].strip(" '")) for k in keyValuePairs]
            keyValues = [k[1].strip() for k in keyValuePairs]
            for i, k in enumerate(keyValues):
                if k == '':
                    keyValues[i] = ''
                else:
                    keyValues[i] = float(k)
            dictOut = dict(zip(keys, keyValues))
        except:
            dictOut = {}
        return dictOut

    def convertString2IntList(self, strInp):
        # Converts a string to an integer list
        try:
            listOut = strInp.strip('[](){}').replace(',', ' ').split()
            if all([k.isdigit() for k in listOut]):
                listOut = [int(k) for k in listOut]
                # Remove repeated items and sort
                acclist = []
                [acclist.append(x) for x in listOut if x not in acclist]
                listOut = sorted(acclist)
            else:
                listOut = None
            if listOut == []:
                listOut = None
        except:
            listOut = None
        return listOut

    def convertDictKeyValues2Strings(self, dict0):
        newDict = {}
        for k in dict0.keys():
            newDict[k] = str(dict0[k])
        return newDict

    def getWestPacificStations(self):
        westPacificStns = {}

        country_filter = \
           ['Australia', 'East Timor', 'Palau', 'Solomon Islands', 'Vanuatu', 'Nauru', 'Niue',
            'Cook Islands', 'Tuvalu', 'Papua New Guinea', 'Marshall Islands', 
            'Micronesia, Federated States of', 'Kiribati', 'Tonga', 'Philippines', 
            'New Caledonia', 'Samoa', 'Fiji']

        for stn in self.stnDict:
            stnlon = self.stnDict[stn][0]
            stnlat = self.stnDict[stn][1]
            #if (stnlon>110) and (stnlon<210) and (stnlat>-25) and (stnlat<10) and (abs(stnlat) > 5):
            #if 1 == 1:
            if (stnlat>-25) and (stnlat<25) and (abs(stnlat) > 5):
                stnCountry = self.stnDict[stn][3]
                stnName = self.stnDict[stn][2]
                if stnCountry in country_filter:
                    if stnCountry in westPacificStns:
                        westPacificStns[stnCountry][stnName] = [stn, stnlon, stnlat]
                    else:
                        westPacificStns[stnCountry] = {}
                        westPacificStns[stnCountry][stnName] = [stn, stnlon, stnlat]
        return westPacificStns       


def MsgDlg(window, string, caption='wxProject', style=wx.YES_NO|wx.CANCEL):
    """Common MessageDialog."""
    dlg = wx.MessageDialog(window, string, caption, style)
    result = dlg.ShowModal()
    dlg.Destroy()
    return result

def getIconData():
    """
    TCRM icon
    License: CC Attribution (http://creativecommons.org/licenses/by/3.0/)
    """
    return zlib.decompress( 
'x\xda\x01\xf2\x03\r\xfc\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00 \x00\
\x00\x00 \x08\x06\x00\x00\x00szz\xf4\x00\x00\x00\x04sBIT\x08\x08\x08\x08|\
\x08d\x88\x00\x00\x03\xa9IDATX\x85\xc5\x97Mh\\U\x14\xc7\x7f\xe7\xbe735\xa36\
\xb14\x934\xa4\xc54\xfdP4j\xa1\x16\xa4\xa1\x9b"(Q\x11\xa9\x01w\xda\x85"\xb5.\
\\\xb8\x10?\xe8F\xd0.,B\xfd\x80\x80T\xdc\x14\x8bY\xd4\x9d\xf1\x8bZ7U(j\x89\
\xa1*MIc\x9aN\x9c\xa4\xd3\xc9\xcc\x9bw\xdfu\xf1&3y\xef\xcdKg:\x81\x1c\xb8\
\x8b\xfbq\xce\xf9\xdd\xff=\xf7>\x9e\xf0\xda1\xc3\x1a\x9aM2Y\xeb9\xe5\xe6\xbc\
\x93\x89\xd5\x00\xa8\x04)k\xcc\xd1W\x1av\xcc\x97]\xee8\xf2\x19X\xb2J\x00J5\
\xe5\x98/\xbb\x90J\x80j\x15 Q\x01\x90\xe6\x00\x00H\xac\x02\xc0\x9d)\xbf\x06\
\x8aJ7\xe5\x98\xb2,\xdaRI\xec\x16\x01d\xe7\xa7\xa7\x8d\xf6<\x04\xb8Qv\xab\
\x13\x06\x98:\xf4t\xb5_p5w\x1f\x1f%Q9*\x01\xda\x126\xb7z\x85\x94\x08Zk\xec\
\xf1?\xfedCG;;z7\x91Y\x97B{\x1e \x91\xc0J\x84-\xe94\x89\x16v,\x80m)\xae\xcd\
\xe7\x19\xbft\xd9\x07@)\xb2\xf3\x0b\x9c\xcd\xcd#"l\xeb\xce\xb0}S\'I\xdb\x0e\
\xb9\x1b6\xa6\x12MKn)\x85\x12av!\xcfo\x93S\xfc\xb7\x90\x07\xabRo"\xd4\xb2\
\x88\xbf\xeb\x89\xe9\x19&\xae\xfc\x0b\xda\xc0\xf0\xbe\xea\xb4B\xd1\x9dn\xc3V\
\xe0\x190\xa6\xa6\x91\x00"B\xa1\xecb\x8c\xc1R\xc2\xc5\xabY.^\xcdr\xedz\x1e\
\xcf\xd5\xb5\xa4V\xb0\xd8\xc3\xdb\xac\xc2\x84\xef\xb7\xa3]F\xc6\xce\xf8\x01\
\x94B*;\x03\xd0e\x97\xa7\x06\xfa\xf9\xea\xe0\x13\x08A\xbf\xfd\xc7O161\x19\
\xabP\xc3wOD\xc0\xb2\xfc&\x821\x06\xedyh\xa7\xcc\xf7\xaf\x0e3z\xf0\xc9Hr\x80\
o^~\x86\x93\xcf\x0f\x81\xe7\xb5\x06\x10g\xef\x0c\xede__\xcf\x8ak\x0e<\xb0\
\x8d\x17\x1e\x19\x80:\x80\xad\x01\x14\x1d\xde~tOd\xb8P\xe7\x9b2\xf2\xec~(\
\x96"\xe3\xf5k\xa0A\xdb\xd5\xdf\x1b\x19\x93\x17\xdf\x85u)(\x951\x1f\xbf\x1e\
\x98\xeb\xe9\xda\xc0T\xeez`\xac%\x05\x1e\xde\x9c\t\xf4\xdf\xff\xee\x17H\xdf\
\xe6\x17j[\x8a\x0f\xcf\x9cg6\xbfXm;;;"1ZR\xc0R\xe1\x9b\x12|\xce\x0f\x9f\xfa\
\x96\xc3_\x8e\xd5\x06DE\xca\xa0%\x80\x9b\x9a\x88\xdfV0\x1b\xa7X\x7f&\xfc\x16\
\x1b\xc0Y\xf4\xaf\xe1\x92\xb9np\x8d\xd6\xc4\xc6\x8b\x05(5\x08\x80\x81b\x11\
\xec\xe5\x00\xa1jw]b\xe3\xb5\x0e\x80\x1f\\/\x03\xd0!\x05\xbc[\x01\x88=\x82\
\x10\x81\xc1\x97\xd7\xbb\x89\x02M\x1f\x81\x13}\x1c\xea\x02\xe0\x81SZ\x19@\
\xbb\xc4\xc6\x8b\x07hP\x01\xc4/Bo\xd9\x074t\x04^H\x81\xbe\xae\x8dt\xdf\xd5\
\xbe\xe4\xcd\x85\xc9+\xcc\xe5o\x84\x00bk\xa0\xce5(\x15A\xd7\x00.O\xcf\x04V\
\x0c\x0f\xee\xe6\xadON\x80m\x83\xe3p\xf6\xf3\x0f\xc8t\xac\xaf\xce\xef>\xf4\
\x06s\xd9l\xc0G\xe1\x94\x88m\xa1\xfc\xe1\xf9\xaf\x7f\xf8)\xb0d{O\x17?\x1e;\
\xc2\xd0C\xf7\xf2\xeb\xc8\xd1@r\x80s\xe7\xceGb\xd8\x94\x1bU\x00_^SS\xc0u53s9\
2\x15\x99\x01\x06\x07\xeea\xf0\xbd7#\xae\xbf\xff}\t\xd0\x84\xf35Q\x84\xc6\
\xa76\xc1\xe7\xb6k\xefc\x98\x0b?\xd7\x8f\xb1\xcc\xee\x7f|\x18nOG\xc6\x15N\
\x91\xd8\x16\xca\x8fS\x88\xae\xb1\x15\xd2\xbf\x8b\xf9Pq-\xd9\xf4l\x16\xe9{\
\x10\x92V\xdd\x1c6\xa5\xc5X\x05N\x8c\x9e\xaevK\x8e\x03\x8e\x03^\x9d\xff\x07\
\x0b\xda\xef\xdbCg\xa6\x93\x97\x9e;\xc0\xd6-\xbd\x8c\xff\xf5\x0f\x1f}q\x92\\\
.\xe7\xff\xc0\xc4\xe4\x116\xefX\xe3\xbf\xe3B~-\xf3\xf3?\x14\xf5\xa0\x97\x88\
\x1c\x93\x19\x00\x00\x00\x00IEND\xaeB`\x82\x17X\xd5\xf9' )

def getIconImage():
    stream = cStringIO.StringIO(getIconData())
    return ImageFromStream(stream)

def getIcon():
    icon = EmptyIcon()
    icon.CopyFromBitmap(BitmapFromImage(getIconImage()))
    return icon


class HelpPopUp ( wx.Dialog ):

    def __init__( self, parent, helpText, boxTitle="Help" ):
        wx.Dialog.__init__  ( self, parent, id = wx.ID_ANY, title = boxTitle, pos = wx.DefaultPosition, size = wx.Size( 391,366 ), style = wx.DEFAULT_DIALOG_STYLE )

        self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )

        bSizer1 = wx.BoxSizer( wx.VERTICAL )

        self.m_panel1 = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.m_panel1.SetFont(parent.guiFont)
        bSizer8 = wx.BoxSizer( wx.VERTICAL )

        bSizer10 = wx.BoxSizer( wx.VERTICAL )

        self.m_richText1 = wx.richtext.RichTextCtrl( self.m_panel1, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.TE_READONLY|wx.RAISED_BORDER|wx.WANTS_CHARS )
        bSizer10.Add( self.m_richText1, 1, wx.EXPAND |wx.ALL, 15 )
        self.m_richText1.BeginAlignment(wx.TEXT_ALIGNMENT_LEFT)

        self.m_richText1.SetFont( parent.guiFont )
        self.m_richText1.WriteText(helpText)

        self.m_button1 = wx.Button( self.m_panel1, wx.ID_ANY, u"Ok", wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer10.Add( self.m_button1, 0, wx.ALIGN_CENTER|wx.ALL, 10 )
        self.m_button1.Bind(wx.EVT_BUTTON, self.okbuttonpressed)
        self.m_button1.SetFocus()

        bSizer8.Add( bSizer10, 4, wx.EXPAND, 5 )

        self.m_panel1.SetSizer( bSizer8 )
        self.m_panel1.Layout()
        bSizer8.Fit( self.m_panel1 )
        bSizer1.Add( self.m_panel1, 1, wx.EXPAND |wx.ALL, 0 )

        self.SetSizer( bSizer1 )
        self.Layout()

    def okbuttonpressed(self, event):
        self.EndModal(wx.ID_OK)


class UnitSelectionBox ( wx.Dialog ):

    def __init__( self, parent, unitName, selectionLabels, initialSelectionStr=""):
        wx.Dialog.__init__  ( self, parent, id = wx.ID_ANY, title = u"TCRM Configuration Editor", pos = wx.DefaultPosition, size = wx.Size( 235,257 ), style = wx.DEFAULT_DIALOG_STYLE )
        self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )

        bSizer1 = wx.BoxSizer( wx.VERTICAL )

        self.m_panel1 = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.m_panel1.SetFont( parent.guiFont )
        bSizer8 = wx.BoxSizer( wx.VERTICAL )

        self.m_staticText1 = wx.StaticText( self.m_panel1, wx.ID_ANY, u"Please select column units:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.m_staticText1.Wrap( -1 )
        bSizer8.Add( self.m_staticText1, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 15 )

        bSizer10 = wx.BoxSizer( wx.VERTICAL )

        self.m_radioBox1 = wx.RadioBox( self.m_panel1, wx.ID_ANY, unitName.rsplit('Units')[0] + ' Units', wx.DefaultPosition, wx.DefaultSize, selectionLabels, 1, wx.RA_SPECIFY_COLS )
        self.m_radioBox1.SetSelection( 0 )
        bSizer10.Add( self.m_radioBox1, 0, wx.ALIGN_CENTER|wx.ALL, 5 )

        self.m_button1 = wx.Button( self.m_panel1, wx.ID_ANY, u"Ok", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.m_button1.Bind(wx.EVT_BUTTON, self.okbuttonpressed)
        bSizer10.Add( self.m_button1, 0, wx.ALIGN_CENTER|wx.ALL, 10 )

        bSizer8.Add( bSizer10, 4, wx.EXPAND, 5 )

        self.m_panel1.SetSizer( bSizer8 )
        self.m_panel1.Layout()
        bSizer8.Fit( self.m_panel1 )
        bSizer1.Add( self.m_panel1, 1, wx.EXPAND |wx.ALL, 0 )

        self.SetSizer( bSizer1 )
        self.Layout()
        if initialSelectionStr is None:
            initialSelectionStr = ""
        # Case-insensitive search for initialSelectionStr in selectionLabels list
        selectionMatch = [k for k in selectionLabels if k.lower() == initialSelectionStr.lower()]
        if len(selectionMatch) == 0:
            initialSelectionStr = ""
        else:
            initialSelectionStr = selectionMatch[0]
        self.m_radioBox1.SetStringSelection(initialSelectionStr)

    def okbuttonpressed(self, event):
        self.EndModal(self.m_radioBox1.GetSelection())


class DateFormatDialog ( wx.Dialog ):

    def __init__( self, parent ):
        wx.Dialog.__init__( self, parent, id = wx.ID_ANY, title = u"TCRM Configuration Editor", pos = wx.DefaultPosition, size = wx.Size( 391,366 ), style = wx.DEFAULT_DIALOG_STYLE )
        self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
        bSizer1 = wx.BoxSizer( wx.VERTICAL )
        self.m_panel1 = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
        self.m_panel1.SetFont( parent.guiFont )
        bSizer8 = wx.BoxSizer( wx.VERTICAL )
        sbSizer2 = wx.StaticBoxSizer( wx.StaticBox( self.m_panel1, wx.ID_ANY, u"Custom Date Format" ), wx.VERTICAL )

        bSizer9 = wx.BoxSizer( wx.HORIZONTAL )
        self.m_staticText1 = wx.StaticText( self.m_panel1, wx.ID_ANY, u"Please specify custom date format:", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.m_staticText1.Wrap( -1 )
        bSizer9.Add( self.m_staticText1, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )
        inputsource = parent.getConfigKeyValue('DataProcess', 'Source')
        if not parent.configDict[inputsource].has_key('DateFormat'):
            parent.setConfigKeyValue(parent.defaultSourceName, 'DateFormat', '%Y-%m-%d %H:%M:%S', valueDataType='str')

        dateFormat = parent.getConfigKeyValue(inputsource, 'DateFormat')

        self.m_textCtrl1 = wx.TextCtrl( self.m_panel1, wx.ID_ANY, dateFormat, wx.DefaultPosition, wx.DefaultSize, 0 )
        bSizer9.Add( self.m_textCtrl1, 1, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5 )

        sbSizer2.Add( bSizer9, 1, wx.EXPAND, 0 )

        self.m_button1 = wx.Button( self.m_panel1, wx.ID_ANY, u"Ok", wx.DefaultPosition, wx.DefaultSize, 0 )
        self.m_button1.Bind(wx.EVT_BUTTON, self.okbuttonpressed)
        sbSizer2.Add( self.m_button1, 0, wx.ALIGN_RIGHT|wx.ALL, 5 )

        bSizer8.Add( sbSizer2, 3, wx.ALL|wx.EXPAND, 15 )
        bSizer10 = wx.BoxSizer( wx.VERTICAL )

        self.m_richText1 = wx.richtext.RichTextCtrl( self.m_panel1, wx.ID_ANY, "", wx.DefaultPosition, wx.DefaultSize, wx.TE_READONLY|wx.RAISED_BORDER|wx.WANTS_CHARS )
        bSizer10.Add( self.m_richText1, 1, wx.EXPAND |wx.ALL, 15 )
        self.m_richText1.BeginAlignment(wx.TEXT_ALIGNMENT_CENTRE);
        self.m_richText1.WriteText(" ")
        self.m_richText1.BeginBold()
        self.m_richText1.WriteText("Help")
        self.m_richText1.EndBold()
        self.m_richText1.EndAlignment();
        self.m_richText1.Newline();
        self.m_richText1.BeginAlignment(wx.TEXT_ALIGNMENT_LEFT);
        self.m_richText1.WriteText("As you have selected a custom date column, you must specify a date format.\n")
        self.m_richText1.WriteText("For example: the format '%Y-%m-%d %H:%M:%S' would accept dates such as '1997-12-30 22:30:00'\n\n")
        self.m_richText1.WriteText("List of date codes:\n")
        self.m_richText1.WriteText(" %Y - Year with century as a decimal number\n")
        self.m_richText1.WriteText(" %y - Year without century as a decimal number [00,99]\n")
        self.m_richText1.WriteText(" %m - Month as a decimal number [01,12]\n")
        self.m_richText1.WriteText(" %d - Day of the month as a decimal number [01,31]\n")
        self.m_richText1.WriteText(" %H - Hour (24-hour clock) as a decimal number [00,23]\n")
        self.m_richText1.WriteText(" %M - Minute as a decimal number [00,59]\n")
        self.m_richText1.WriteText(" %S - Second as a decimal number [00,59]\n")
        bSizer8.Add(bSizer10, 4, wx.EXPAND, 10)
        self.m_panel1.SetSizer(bSizer8)
        self.m_panel1.Layout()
        bSizer8.Fit(self.m_panel1)
        bSizer1.Add(self.m_panel1, 1, wx.EXPAND |wx.ALL, 0)
        self.parentpointer = parent
        self.SetSizer( bSizer1 )
        self.Layout()

    def okbuttonpressed(self, event):
        inputsource = self.parentpointer.getConfigKeyValue('DataProcess', 'Source')
        self.parentpointer.setConfigKeyValue(inputsource, 'DateFormat', self.m_textCtrl1.GetValue().strip("\" '"))
        self.EndModal(wx.ID_OK)


if __name__ == "__main__":
	app = wx.PySimpleApp(0)
	wx.InitAllImageHandlers()
	toplevelwin = MainFrame(None)
	app.SetTopWindow(toplevelwin)
	toplevelwin.Show()
	app.MainLoop()
