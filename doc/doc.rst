======================
PyZgoubi Documentation
======================

The latest version of PyZgoubi and this document can be found at http://www.hep.man.ac.uk/u/sam/pyzgoubi/

.. contents::


Introduction
------------

http://sourceforge.net/projects/zgoubi Zgoubi is a particle tracking code maintained by François Méot. PyZgoubi is an interface to Zgoubi writing in python. It aims to make input files that are easy to read and can contain calculations, loops, and any other python program feature.

A basic knowledge of http://python.org/ Python is useful to use PyZgoubi. Python is a general purpose, high level, interpreted programing language. I recommend the http://www.python.org/doc/2.5.2/tut/tut.html Official Python Tutorial.

Quick Start
-----------

Put the following in a text file named quickstart.py::

    make_line('line')

    # create and OBJET2 and add 1 electron to it
    ob = OBJET2()
    ob.set(BORO=ke_to_rigidity(10e6, ELECTRON_MASS))
    ob.add(Y=0, T=0.1, D=1)
    add(ob)

    add(FAISCNL(FNAME='zgoubi.fai')) # record point
    d1 = DRIFT(XL=50)
    q1 = QUADRUPO(XL=20, B_0=5, R_0=10, XPAS=(10,20,10))
    d2 = DRIFT(XL=50)
    q2 = QUADRUPO(XL=20, B_0=-5, R_0=10, XPAS=(10,20,10))
    add(d1, q1, d2, q2)
    add(FAISCNL(FNAME='zgoubi.fai')) # record point
    add(END())

    print output() # file fed to zgoubi (zgoubi.dat)
    run() # run the line in zgoubi
    print res() # show the zgoubi.res file
    print get_all('fai') # so the data recorded with FAISCNL

The file can also be found in examples/quickstart.py

Run it with::

    pyzgoubi quickstart.py

You should see the zgoubi data file, the output from zgoubi, the res file from zgoubi, and the coordinates at the 2 FAISCNL points displayed to the terminal.

For more advance use it may be better store the Line with a name, this allows you to work with multiple Lines simultaneously. The following shows the same file, but with a named line::

    myline = Line('line')

    # create and OBJET2 and add 1 electron to it
    ob = OBJET2()
    ob.set(BORO=ke_to_rigidity(10e6, ELECTRON_MASS))
    ob.add(Y=0, T=0.1, D=1)
    myline.add(ob)

    myline.add(FAISCNL(FNAME='zgoubi.fai')) # record point
    d1 = DRIFT(XL=50)
    q1 = QUADRUPO(XL=20, B_0=5, R_0=10, XPAS=(10,20,10))
    d2 = DRIFT(XL=50)
    q2 = QUADRUPO(XL=20, B_0=-5, R_0=10, XPAS=(10,20,10))
    myline.add(d1, q1, d2, q2)
    myline.add(FAISCNL(FNAME='zgoubi.fai')) # record point
    myline.add(END())

    print myline.output() # file fed to zgoubi (zgoubi.dat)
    myresults = myline.run() # run the line in zgoubi
    print myresults.res() # show the zgoubi.res file
    print myresults.get_all('fai') # so the data recorded with FAISCNL 

Generating zgoubi.dat files
---------------------------

The simplest level of Pyzgoubi is to use it to generate a zgoubi.dat file. You can then run Zgoubi with the file to get your results.

PyZgoubi works with Lines, which are made up of Elements. The Elements are mostly the same as those used in Zgoubi. Most Elements correspond to magnets and other beam line elements of a physical accelerator, e.g. DIPOLE, MUTIPOL and CAVITE. Some are used to define the beam e.g. OBJET2. Some are to control the execution of Zgoubi, e.g. FAISCNL, MATRIX and END. A number of these Elements can be added to a Line. This Line then has all the information needed to make a zgoubi.dat file. 

A Line can be created as follows::

    emma = Line('emma simulation')

emma is the name used to refer to the line in the program. 'emma simulation' is the text that will be put at the start of the zgoubi.dat file.

To create an Element use::

    q1 = QUADRUPO('defoc', XL=20, R_0=2, B_0=2, XPAS=0.1)

This can now be added to the Line ::

    emma.add(q1)

These last 2 steps can be contracted::

    emma.add(QUADRUPO('defoc', XL=20, R_0=2, B_0=2, XPAS=0.1))

however this means that there is no reference to the element that could be used for modifying its parameters later, e.g.::

    q1.set(B_0=3)

The line can now be used to output a zgoubi.dat using the output() function::

    print emma.output()

All these instructions can be put in a text file and run using the command pyzgoubi. (pyzgoubi is an alias set up by the installer). Below is a section of emma.py from the examples::

        emma = Line('emma')
        xpas = (20,20,20)

        cells = 42
        angle = 360/cells
        d_offset = -34.048 * mm
        f_offset = -7.514 * mm

        #lengths
        ld = 210 * mm
        sd = 50 * mm

        fq = 58.782 * mm
        dq = 75.699 * mm

        # quad radius
        fr = 37 * mm
        dr = 53 * mm

        fb = -6.695 * fr * T
        db = 4.704 * dr * T

        ob = OBJET2()
        emma.add(ob)

        emma.add(ELECTRON())

        emma.add(DRIFT('ld', XL=ld*cm_/2))
        emma.add(CHANGREF(ALE=angle))

        emma.add(CHANGREF(YCE=d_offset*cm_))
        emma.add(QUADRUPO('defoc', XL=dq*cm_, R_0=dr*cm_, B_0=db*kgauss_, XPAS=xpas))

        emma.add(CHANGREF(YCE=-d_offset*cm_))

        emma.add(DRIFT('sd', XL=sd*cm_))

        emma.add(CHANGREF(YCE=f_offset*cm_))
        emma.add(QUADRUPO('foc', XL=fq*cm_, R_0=fr*cm_, B_0=fb*kgauss_, XPAS=xpas))
        emma.add(CHANGREF(YCE=-f_offset*cm_))

        emma.add(DRIFT('ld', XL=ld*cm_/2))

        emma.add(FAISCNL(FNAME='zgoubi.fai'))

        emma.add(REBELOTE(K=99, NPASS=10))

        emma.add(END())

        rigidity = ke_to_rigidity(10e6, 0.51099892e6)
        ob.set(BORO=-rigidity)
        ob.add(Y=0, T=0, D=1)

        print emma.output()

This can be run with the command::

    pyzgoubi emma.py


It will build a Line and print the zgoubi.dat input to the screen. The '.py' extension is not necessary, but will cause your text editor to use python syntax highlighting.



Elements
--------

There are two ways Elements can be defined in pyzgoubi. Most Elements are simple, they have a static list of parameters. Some have some extra complexity, for example different parameters depending on options, sections repeated N times. These elements can be defined using a simple syntax, which is then converted into python code. More complex elements must be written in python.


The simple elements are defined in defs/simple_elements.defs. For each element there is a number of lines of text, delimited by blank lines. Comments can be put after a '#' character. The first line gives the name of the class, this is the name you use in the input file. The second line gives the name used by Zgoubi, this must match the Zgoubi manual. Then follows a line for each line of output in the zgoubi.dat file; first the names of the parameters, then a ':', then the types. For example::


    BEND
    BEND
    IL : I
    XL, Sk, B1 : 3E
    X_E, LAM_E, W_E : 3E
    N, C_0, C_1, C_2, C_3, C_4, C_5 : I,6E
    X_S, LAM_S, W_S : 3E
    NS, CS_0, CS_1, CS_2, CS_3, CS_4, CS_5 : I,6E
    XPAS: X
    KPOS, XCE, YCE, ALE : I,3E

The types can be:

- I : integer
- E : real (floating point)
- Ax : string with up to x characters
- X : special type for XPAS. Can be integer, or group of 3 integers e.g. (10,20,10)

The type can be followed by a number for several parameters of the same type.

If the parameters used vary depending on the value of another option the following syntax can be used::

    CAVITE
    CAVITE
    IOPT : I
    !IOPT==0
    X, X : 2E
    !IOPT==1
    L, h : 2E
    V, X : 2E
    !IOPT==2
    L, h : 2E
    V, sig_s : 2E
    !IOPT==3
    X,X : 2E
    V, sig_s : 2E

Here the value of IOPT switches the element to output different parameters. (See the zgoubi manual's description of CAVITE for more info).


Elements with a looped section can be defined as follows::

    FFAG
    FFAG
    IL : I
    N, AT, RM: I,2E
    !N*{
    ACN, DELTA_RM, BZ_0, K: 4E
    G0_E, KAPPA_E: 2E
    NCE, CE_0, CE_1, CE_2, CE_3, CE_4, CE_5, SHIFT_E: I,7E
    OMEGA_E, THETA_E, R1_E, U1_E, U2_E, R2_E: 6E
    G0_S, KAPPA_S: 2E
    NCS, CS_0, CS_1, CS_2, CS_3, CS_4, CS_5, SHIFT_S: I,7E
    OMEGA_S, THETA_S, R1_S, U1_S, U2_S, R2_S: 6E
    G0_L, KAPPA_L: 2E
    NCL, CL_0, CL_1, CL_2, CL_3, CL_4, CL_5, SHIFT_L: I,7E
    OMEGA_L, THETA_L, R1_L, U1_L, U2_L, R2_L: 6E
    !}
    KIRD, RESOL: 2I
    XPAS: E
    KPOS, RE, TE, RS, TS: I,4E

Here the section between the braces is repeated. These elements are used slightly differently to simpler elements. Then non looping section is defined normally.::


    triplet = FFAG('triplet', IL=0, AT=10 ... )

Then the looped part can be added::

    triplet.add(ACN = 6, BZ_0 = 0.5 ...)
    triplet.add(ACN = 4, BZ_0 = -0.5 ...)
    triplet.add(ACN = 6, BZ_0 = 0.5 ...)

N gets automatically set. All the looped parts can be removed using::

    triplet.clear()

When pyzgoubi runs it searches the defs folder for files ending in .defs. Additional files can be added to the extra_defs_files list in zgoubi_settings.py. If any of these files have been modified then they are reread and the defs.py is regenerated.


The Elements that cannot be defined in this way must be put into the static_defs.py file. They must be classes that have an output() method, which generates the code needed for the zgoubi.dat file.


There is also a FAKE_ELEM element. This allows you to put arbitrary text into the zgoubi.dat file. It is useful for using an Zgoubi element that pyzgoubi does not have a definition for. For example::

    change_txt = """'CHANGREF'
    5.0 0 10.0
    """
    change = FAKE_ELEM(change_txt)
    line.add(change)



Available Elements
""""""""""""""""""

- BEND
- CAVITE
- CHANGREF
- DRIFT
- ELECTRON
- END
- FAISCEAU
- FAISCNL
- FAISTORE
- FAKE_ELEM
- FFAG
- IMMORTAL_MUON
- MARKER
- MATRIX
- MULTIPOL
- MUON
- OBJET1
- OBJET2
- OBJET5
- PARTICUL
- PROTON
- QUADRUPO
- REBELOTE
- TOSCA

To find the full list of elements available in the current version run::

    pyzgoubi help elements

To find the names of the parameters available for an element use::

    pyzgoubi help element_name

e.g.::

    pyzgoubi help MULTIPOL

Use this in combination with the Zgoubi manual.

Running Zgoubi
--------------

Once a line has been created and had the needed elements added it can be run. PyZgoubi will take care of creating a temporary directory, creating the zgoubi.dat file and running Zgoubi. This is done to prevent zgoubi from overwriting any existing files. If you wish to keep any of the output files you must use the commands to copy these to where you want them.


The following example shows how to run a line::


    #create line
    emma = Line('emma')

    #add elements
    emma.add( ...  )
    ...

    #run line
    emma.run()

    #save output
    emma.save_res("emma.res")
    emma.save_plt("emma.plt")

Note that you will need to make sure your line will actually create plt or fai files, otherwise you will receive a file not found error. See the Zgoubi manual for more information.

The run command can take several options. If you want to inspect the directory where zgoubi is run, or to use zpop, then set xterm=True. If you want to change the directory that zgoubi is run in you can use the tmp_prefix option. It is best to make sure this is a local disk (i.e. not a network/remote disk). The default directory can be set in the zgoubi_settings.py file.::

    emma.run(xterm=True)
    emma.run(tmp_prefix = '/var/tmp/sam/')
    emma.run(xterm=True, tmp_prefix = '/home/sam/tmp/')

If you want to do analysis of the simulation you can use the Results object that is returned by the run() function.::

    res = emma.run()

See the Results Object chapter for more info.

Each time a line is run a temporary director is created. These are normally automatically cleared up when PyZgoubi finishes (also the /tmp directory is usually emptied when a computer shuts down). However if you are making repeated calls to run(), then you may want to manually clear away these files. This can be done with the clean() function. This removes all the temporary directories the currently running PyZgoubi has made for the line. Don't clean the line until you have finished working with its output files.::

    emma.run()
    emma.clean()


Results Object
--------------

When you run a line it creates a Results object, that can be used to get information about the paths of the particles.::

    res = pamela.run()

get_all() and get_track(), let you get lists of the particle coordinates. They each need to be told if they should read the plt (points within the magnetic elements) or fai (beam at FAISCNL element). get_all() returns a list of dictionaries, containing all the coordinates and information. get_track() returns a list of lists of just the requested coordinates.::

    print res.get_all('plt')
    print res.get_all('fai')
    print res.get_track('fai', ['Y','T'])


Loops
-----

For making complex lines it can be useful to use python features such as loops, e.g. to put 5 identical FODO cells you could use ::

    line = Line("example")

    line.add( ... )
    line.add( ... )
    
    for x in range(5):
        line.add(QUADRUPO( ... ))
        line.add(DRIFT( ... ))
        line.add(QUADRUPO( ... ))
        line.add(DRIFT( ... ))

    line.add( ... )

Note that the for loop block lasts as long as the code is indented.

If you want to make one iteration different then you can do a test based on the x ::

    for x in range(5):
        line.add(QUADRUPO( ... ))
        line.add(DRIFT( ... ))
        line.add(QUADRUPO( ... ))
        line.add(DRIFT( ... ))
        if (x == 2): #note x counts from zero
            line.add(DIPOLE( ... )

Command line arguments
----------------------
There is a convenience function for using command line arguments. For example near the start of the code put::

    number_of_laps = int(get_cmd_param('laps', 10))

Then when creating the line use::

    l.add(REBELOTE(NPASS=number_of_laps-1, KWRIT=1, K=99))

to access the variable. When you run the your simulation use the command line argument::

    pyzgoubi sim.py laps=50

I the second parameter to get_cmd_param() is the default. If you don't give a value as an argument then the default is used. If you don't give a default eg::

    particle_energy = float(get_cmd_param('energy'))

then you will receive and error if you don't provide a value as an argument.

Tips
----
Python hints
""""""""""""
If there is a '#' character on a line, everything after it is treated as a comment.

Python uses whitespace to delimit blocks (instead of braces '{' and '}' in C/C++). The PyZgoubi code uses tabs, so it is best to use tabs in your input files. If you get an 'IndentationError' check that you have not mixed spaces and tabs, or accidentally started a line with a space/tab.

Python identifiers (variable, function, object names etc) are case sensitive, they must start with a letter and only contain letters, numbers and underscores.



