#import subprocess
#import sys

#def install(package):
    #subprocess.check_call([sys.executable, "-m", "pip", "install", package])

#install('tkinter')
#install('pillow')
import threading
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkmacosx import Button
import subprocess as sub
import time
import os
import sys
#import pmw



splash_root = Tk()

img = Image("photo", file="Mini.png")
splash_root.call('wm','iconphoto', splash_root._w, img)
img2 = PhotoImage(file="Splash.png")

#bg=PhotoImage(file="Ilustración_sin_título copia.png")

screen_width = splash_root.winfo_screenwidth()
screen_height = splash_root.winfo_screenheight()
app_width = 480
app_height = 480
x = (screen_width / 2) - (app_width / 2)
y = (screen_height / 2) - (app_height / 2)

x=round(x)
y=round(y)

canvas = Canvas(splash_root, width=400, height=400, bg='white', highlightthickness=0)
canvas.pack()
canvas.master.overrideredirect(True)
canvas.master.geometry('+' + str(x) + '+' + str(y) + '')
canvas.master.wm_attributes("-transparent", "true")
canvas.master.wm_attributes("-topmost", True)
canvas.master.lift()
mylabel=Label(splash_root, image=img2)
mylabel.place(x=0, y=0, relwidth=1, relheight=1)

w = 0  # index number

#To not get errors
SS_model1 = ''
S_model1 = ''
SStadistics = ''

def main_window():
    splash_root.destroy()

    root = Tk()
    root.title('ProtModel')
    root.geometry("500x480")
    root.configure(bg='gray90')
    root.resizable(False, False)
    root.bind("<Button-1>", lambda event: event.widget.focus_set())
    root.bind("<Button-1>", root.focus_set())

    s = ttk.Style()
    s.theme_use('clam')
    s.theme_settings('clam', {
        "TCombobox": {
            "configure": {"padding": 1},
            "map": {
                "background": [("active", "lightblue"),
                               ("!disabled", "lightblue")],
                "fieldbackground": [("!disabled", "white")],
                "foreground": [("focus", "black"),
                               ("!disabled", "black")],
                "arrowcolor": [("focus", "darkblue"),
                               ("!disabled", "darkblue")],
                "bordercolor": [("focus", "darkblue"),
                               ("!disabled", "darkblue")],
                "darkcolor": [("focus", "darkblue"),
                                ("!disabled", "darkblue")],
                "focusfill": [("focus", "darkblue"),
                                ("!disabled", "darkblue")],
                "selectforeground": [("focus", "black"),
                              ("!disabled", "black")],
                "selectbackground": [("focus", "white"),
                                     ("!disabled", "white")],
            }
        }
    })

    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    app_width = 500
    app_height = 480

    x = (screen_width / 2) - (app_width / 2)
    y = (screen_height / 2) - (app_height / 2)

    root.geometry(f'{app_width}x{app_height}+{int(x)}+{int(y)}')

    global is_on
    is_on = True

    main_frame = Frame(root, bg='gray90')
    main_frame.pack(fill=BOTH, expand=1)

    my_canvas = Canvas(main_frame, bg='gray90')
    my_canvas.pack(side=LEFT, fill=BOTH, expand=1)

    my_scrollbary = ttk.Scrollbar(main_frame, orient=VERTICAL, command=my_canvas.yview)
    my_scrollbary.pack(side=RIGHT, fill=Y)
    # my_scrollbarx = ttk.Scrollbar(main_frame, orient=HORIZONTAL, command=my_canvas.xview)
    # my_scrollbarx.pack(side=BOTTOM, fill=X)

    my_canvas.configure(yscrollcommand=my_scrollbary.set)
    # my_canvas.configure(xscrollcommand=my_scrollbarx.set)
    my_canvas.bind('<Configure>', lambda e: my_canvas.configure(scrollregion=my_canvas.bbox("all")))

    second_frame = Frame(my_canvas, bg='gray90')
    my_canvas.create_window((0, 0), window=second_frame, anchor="nw")

    def on_vertical(event):
        my_canvas.yview_scroll(-1 * event.delta, 'units')

    #def on_horizontal(event):
        #my_canvas.xview_scroll(-1 * event.delta, 'units')

    second_frame.bind('<MouseWheel>', on_vertical)
    #my_canvas.bind_all('<Shift-MouseWheel>', on_horizontal)
    combo=ttk.Combobox()
    combo.unbind_class("TCombobox", "<MouseWheel>")
    
    gaps = [
        "Ignored",
        "New State"
    ]

    S_Simulations = [
        "No",
        "Yes"
    ]

    S_Running = [
        "No",
        "Yes"
    ]

    haploid = [
        "",
        "Haploid",
        "Diploid"
    ]

    distribution = [
        "",
        "fix",
        "uniform",
        "norm",
        "exp",
        "gamma",
        "beta",
        "dirichlet"
    ]

    amino_distribution = [
        "Default",
        "fix",
        "dirichlet"
    ]

    truncated = [
        "",
        "t "
    ]

    S_model = [
        "Blosum62",
        "CpRev",
        "Dayhoff",
        "DayhoffDCMUT",
        "HIVb",
        "HIVw",
        "JTT",
        "JonesDCMUT",
        "LG",
        "Mtart",
        "Mtrev24",
        "RtRev",
        "VT",
        "WAG"
        #"UserEAAM"
    ]

    ABCMethod = [
        "",
        "rejection",
        "mnlogistic",
        "neuralnet"
    ]

    MultiPages = [
        "",
        "No",
        "Yes"
    ]

    GrowthRate = [
        " ",
        "Exponential Growth Rate",
        "Demographic periods"
    ]

    MigrationModel = [
        "",
        "0",
        "1",
        "2",
        "3"
    ]

    SubsModSelec = list()

    def popup(tit, messag):
        pop = Toplevel()
        pop.title(tit)
        #pop.geometry("350x125")
        pop.config(bg="white")
        pop.geometry(f"350x125+{root.winfo_x()+70}+{root.winfo_y()+200}")
        pop.resizable(False, False)

        pop_ico = Label(pop, text="❌", font=('Arial', 40), fg='red', bg='white')
        #pop_ic2 = Label(pop, text="⚠", font=('Arial', 40), fg='red', bg='gray90')
        #pop_ic3 = Label(pop, text="🚨", font=('Arial', 40), fg='red', bg='gray90')
        pop_ico.grid(row=1, column=1, pady=10, padx=10)
        pop_label1 = Label(pop, text = messag, font=('Arial', 14), bg='white')
        pop_label1.grid(row=1, column=2, padx=10)
        Boton_pop = Button(pop, text="OK", command=pop.destroy, bg="gray90")
        Boton_pop.grid(row=2, column=2, pady=10)
        
    def popup2(tit, messag):
        pop = Toplevel()
        pop.title(tit)
        #pop.geometry("350x125")
        pop.config(bg="white")
        pop.geometry(f"350x125+{root.winfo_x()+70}+{root.winfo_y()+200}")

        pop_ico = Label(pop, text="🥳", font=('Arial', 40), fg='red', bg='white')
        #pop_ic2 = Label(pop, text="⚠", font=('Arial', 40), fg='red', bg='gray90')
        #pop_ic3 = Label(pop, text="🚨", font=('Arial', 40), fg='red', bg='gray90')
        pop_ico.grid(row=1, column=1, pady=10, padx=10)
        pop_label1 = Label(pop, text = messag, font=('Arial', 14), bg='white')
        pop_label1.grid(row=1, column=2, padx=10)
        Boton_pop = Button(pop, text="OK", command=pop.destroy, bg="gray90")
        Boton_pop.grid(row=2, column=2, pady=10)

    def updateScrollRegion():
        my_canvas.update_idletasks()
        my_canvas.config(scrollregion=second_frame.bbox())

    def openSS():

        newWindows = Toplevel(second_frame)
        newWindows.title("Summary stadistics")
        #newWindows.geometry("285x315")
        newWindows.geometry(f"285x160+{root.winfo_x()+110}+{root.winfo_y()}")
        newWindows.config(background="gray90")
        newWindows.resizable(False, False)
        newWindows_frame = ttk.Frame(newWindows, padding=(6, 3, 12, 12))
        newWindows_frame.grid(column=0, row=0, sticky=(N, S, E, W))
        #s.configure('TFrame', background='gray90')
        # valores
        newWindows_valores = StringVar()
        newWindows_valores.set("DGREM_mean DGREM_sd SegSites Grantham_mean_Position Grantham_sd_Position Grantham_sk_Position Grantham_ku_Position")
        # crear la lista
        newWindows_lstbox = Listbox(newWindows_frame, listvariable=newWindows_valores, selectmode=MULTIPLE, width=30, height=7)
        newWindows_lstbox.select_set(first=0, last=7)
        newWindows_lstbox.grid(column=0, row=0, columnspan=2)
        # frame.bind('<MouseWheel>', no_op)

        # selecionar
        def select():
            global SStadistics
            #SStadistics = ''
            SumStadistics = []
            seleccion = newWindows_lstbox.curselection()
            for i in seleccion:
                entrada = newWindows_lstbox.index(i) +1
                SumStadistics.append(entrada)
            SStadistics = str(SumStadistics)[1:-1]
            SStadistics = SStadistics.replace("'", "")
            SStadistics = SStadistics.replace(",", "")


            # newWindow.destroy()

        Select = Button(newWindows_frame, text="Select", command=select, bg='white', borderless=1, focuscolor='')
        Select.grid(column=0, row=1, padx=85, pady=5)

    def openSM():

        # Toplevel object which will be treated as a new window
        newWindow = Toplevel(second_frame)
        # sets the title of the Toplevel widget
        newWindow.title("Substitution models")
        # sets the geometry of toplevel
        #newWindow.geometry("285x315")
        newWindow.geometry(f"285x375+{root.winfo_x()+110}+{root.winfo_y()}")
        newWindow.config(background="gray90")
        newWindow.resizable(False, False)
        # create a frame
        newWindow_frame = ttk.Frame(newWindow, padding=(6, 3, 12, 12))
        newWindow_frame.grid(column=0, row=0, sticky=(N, S, E, W))
        s.configure('TFrame', background='gray90') #Para mismo color s=estilo
        # valores
        newWindow_valores = StringVar()
        newWindow_valores.set("Blosum62 CpRev Dayhoff DayhoffDCMUT HIVb HIVw JTT JonesDCMUT LG Mtart Mtrev24 RtRev VT WAG UserEAAM")
        # crear la lista
        lab9999 = Label(newWindow_frame, text="Empirical substitution models", bg="gray90")
        lab9999.grid(column=0, row=0, sticky="SW")
        newWindow_lstbox = Listbox(newWindow_frame, listvariable=newWindow_valores, selectmode=MULTIPLE, width=30, height=10)
        newWindow_lstbox.select_set(first=0, last=0)
        newWindow_lstbox.grid(column=0, row=1,sticky="SW")
        #frame.bind('<MouseWheel>', no_op)
        lab99=Label(newWindow_frame, text="Structural substitution models", bg="gray90")
        lab99.grid(column=0, row=2, sticky="SW")
        #v1 = IntVar()

        #def v1_var(value):
            #if v1.get() == 0:
                #Struc_Model_F.configure(state="disabled")
                #Struc_Model_N.configure(state="disabled")
                #v2.set(0)
                #v3.set(0)
            #else:
                #Struc_Model_F.configure(state="normal")
                #Struc_Model_N.configure(state="normal")


        #Struc_Model_Yes = Radiobutton(newWindow_frame, text = "Yes", variable = v1, bg="gray90", value =1, command=lambda: v1_var(v1.get()))
        #Struc_Model_No = Radiobutton(newWindow_frame, text = "No", variable=v1, bg="gray90", value =0, command=lambda: v1_var(v1.get()))
        #Struc_Model_Yes.grid(column=0, row=3, sticky="SW")
        #Struc_Model_No.grid(column=1, row=3, sticky="SW")
        v2 = IntVar()
        v3 = IntVar()
        Struc_Model_F = Checkbutton(newWindow_frame, text = "Fitness", variable = v2, bg="gray90")
        Struc_Model_N = Checkbutton(newWindow_frame, text = "Neutral", variable = v3, bg="gray90")
        Struc_Model_F.grid(column=0, row=3, sticky="SW")
        Struc_Model_N.grid(column=0, row=3, sticky="S")
        Struc_Model_F.select()
        Struc_Model_N.select()

        Chainlab=Label(newWindow_frame, text="Select PDB chain", bg="gray90")
        Chainlab.grid(column=0, row=4, sticky="SW")
        ChainEntry=Entry(newWindow_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
        ChainEntry.config(width=3)
        ChainEntry.grid(column=0, row=4, sticky="S")

        v4 = IntVar()
        def activate_button():
            if v4.get() == 1:
                SelectGMRCA.config(state=NORMAL)
            else:
                SelectGMRCA.config(state=DISABLED)
        GMRCAlab=Label(newWindow_frame, text="Use GMRCA sequence?", bg="gray90")
        GMRCAlab.grid(column=0, row=5, sticky="SW")
        GMRCA_Yes = Radiobutton(newWindow_frame, text = "Yes", variable=v4, bg="gray90", value =1, command=lambda: activate_button())#, command=lambda: v1_var(v1.get()))
        GMRCA_No = Radiobutton(newWindow_frame, text = "No", variable=v4, bg="gray90", value =0, command=lambda: activate_button())#, command=lambda: v1_var(v1.get()))
        GMRCA_Yes.grid(column=0, row=6, sticky="SW")
        GMRCA_No.grid(column=0, row=6, sticky="S")
        global GMRCA1
        GMRCA1 = ''

        def GMRCA_1():
            GMRCA = filedialog.askopenfilename()
            GMRCA1 = GMRCA.split('/')[-1]


        SelectGMRCA = Button(newWindow_frame, width=70, text="GMRCA", command=GMRCA_1, bg='white', borderless=1, focuscolor='', state=DISABLED)
        SelectGMRCA.grid(column=0, row=6, sticky="SE", padx=20)


        saveandcontinue=Label(newWindow_frame, text="Save changes and select the template?", bg="gray90")
        saveandcontinue.grid(column=0, row=7, sticky="SW")

        # selecionar
        def select():
            global pdb_file1
            global chain1
            if v2.get() or v3.get() == 1:
                pdb_file = filedialog.askopenfilename()
                pdb_file1 = pdb_file.split('/')[-1]
                pdb_filename = Label(newWindow_frame, text=pdb_file1, bg="gray90")
                pdb_filename.grid(column=0, row=8, sticky="SW")
            else:
                popup("ERROR", "If you select a template you must select at least one structutal substitution model (Fitness, Neutral)")

            if ChainEntry.get() != '':
                chain1 = ChainEntry.get()
            else:
                popup("ERROR", "You must select the PDB chain")

            global SS_model1
            #SS_model1 = ''
            global SS_model
            SS_model = []
            #global SS_model
            if v2.get() == 1:
                SS_model.append('Fitness')
            if v3.get() == 1:
                SS_model.append('Neutral')
            SS_model1 = str(SS_model)[1:-1]
            SS_model1 = SS_model1.replace("'", "")
            SS_model1 = SS_model1.replace(",", "")


            global S_model1
            #S_model1 = ''
            S_models = []
            seleccion = newWindow_lstbox.curselection()
            for i in seleccion:
                entrada = newWindow_lstbox.get(i)
                S_models.append(entrada)
            S_model1 = str(S_models)[1:-1]
            S_model1 = S_model1.replace("'", "")
            S_model1 = S_model1.replace(",", "")



            #newWindow.destroy()

        Select = Button(newWindow_frame, text="Select", command=select, bg='white', borderless=1, focuscolor='')
        Select.grid(column=0, row=8, sticky="SE", padx=20)

    #global alignment_directory
    def Name_file():
        global alignment_directory
        global alignment_directory2
        global filename
        filename = Label(second_frame, text=" ", bg="gray90")
        alignment_directory = filedialog.askopenfilename()
        if alignment_directory == '':
            alignment_directory2 = ''
        else:
            alignment_directory2 = alignment_directory.split('/')[-1]
            os.chdir(alignment_directory[:-len(alignment_directory2)])
            #filename = Label(second_frame, text=" ", bg="gray90")
        if alignment_directory2 != '':
            lab3_2.config(text=alignment_directory2)
            filename.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            filename.grid(row=4, column=1, sticky="W", padx=5)
        else:
            filename.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            filename.grid(row=4, column=1, sticky="W", padx=5,)
            lab3_2.config(text='')
            popup("ERROR","You must upload an alignment file")

    #def pick_reco(e):
        #if recombination_menu.get() == "fix":
            #recombination_menu_1.delete(first=0,last=1000)
            #recombination_menu_2.current(0)
            #recombination_menu_2.configure(state="disabled")
            #recombination_menu_3.delete(first=0,last=1000)
        #if recombination_menu.get() == "uniform":
            #recombination_menu_1.delete(first=0,last=1000)
            #recombination_menu_2.current(0)
            #recombination_menu_2.configure(state="disabled")
            #recombination_menu_3.delete(first=0,last=1000)
        #if recombination_menu.get() == "norm":
            #recombination_menu_1.delete(first=0,last=1000)
            #recombination_menu_2.current(0)
            #recombination_menu_2.configure(state="enabled")
            #recombination_menu_3.delete(first=0,last=1000)
        #if recombination_menu.get() == "exp":
            #recombination_menu_1.delete(first=0,last=1000)
            #recombination_menu_2.current(0)
            #recombination_menu_2.configure(state="enabled")
            #recombination_menu_3.delete(first=0,last=1000)
        #if recombination_menu.get() == "gamma":
            #recombination_menu_1.delete(first=0,last=1000)
            #recombination_menu_2.current(0)
            #recombination_menu_2.configure(state="enabled")
            #recombination_menu_3.delete(first=0,last=1000)
        #if recombination_menu.get() == "beta":
            #recombination_menu_1.delete(first=0,last=1000)
            #recombination_menu_2.current(0)
            #recombination_menu_2.configure(state="enabled")
            #recombination_menu_3.delete(first=0,last=1000)
        #if recombination_menu.get() == "dirichlet":
            #recombination_menu_1.delete(first=0,last=1000)
            #recombination_menu_2.current(0)
            #recombination_menu_2.configure(state="disabled")
            #recombination_menu_3.delete(first=0,last=1000)
            #recombination_menu_3.configure(state="disabled")
    
    def pick_SR(e):
        if SR_menu.get() == "fix":
            SR_menu_1.delete(first=0,last=1000)
            SR_menu_2.current(0)
            SR_menu_2.configure(state="disabled")
            SR_menu_3.delete(first=0,last=1000)
            SR_menu_3.configure(state="disabled")
        if SR_menu.get() == "uniform":
            SR_menu_1.delete(first=0,last=1000)
            SR_menu_2.current(0)
            SR_menu_2.configure(state="disabled")
            SR_menu_3.delete(first=0,last=1000)
            SR_menu_3.configure(state="disabled")
        if SR_menu.get() == "norm":
            SR_menu_1.delete(first=0,last=1000)
            SR_menu_2.current(0)
            SR_menu_2.configure(state="enabled")
            SR_menu_3.delete(first=0,last=1000)
        if SR_menu.get() == "exp":
            SR_menu_1.delete(first=0,last=1000)
            SR_menu_2.current(0)
            SR_menu_2.configure(state="enabled")
            SR_menu_3.delete(first=0,last=1000)
        if SR_menu.get() == "gamma":
            SR_menu_1.delete(first=0,last=1000)
            SR_menu_2.current(0)
            SR_menu_2.configure(state="enabled")
            SR_menu_3.delete(first=0,last=1000)
        if SR_menu.get() == "beta":
            SR_menu_1.delete(first=0,last=1000)
            SR_menu_2.current(0)
            SR_menu_2.configure(state="enabled")
            SR_menu_3.delete(first=0,last=1000)
        if SR_menu.get() == "dirichlet":
            SR_menu_1.delete(first=0,last=1000)
            SR_menu_2.current(0)
            SR_menu_2.configure(state="disabled")
            SR_menu_3.delete(first=0,last=1000)
            SR_menu_3.configure(state="disabled")

    A_frequenciesModel = 'fix '
    def pick_amino(e):
        if A_frequencies_menu.get() == "Default":
            A_frequencies.configure(state="normal")
            A_frequencies.delete(first=0, last=1000)
            A_frequencies.insert(0, "0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05")
            A_frequencies.configure(state="disabled")
            A_frequenciesModel = 'fix '
        elif A_frequencies_menu.get() == "dirichlet":
            A_frequencies.configure(state="normal")
            A_frequencies.delete(first=0, last=1000)
            A_frequencies.insert(0, "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1")
            A_frequencies.configure(state="disabled")
            A_frequenciesModel = 'dirichlet '
        else:
            A_frequencies.configure(state="normal")
            A_frequencies.delete(first=0, last=1000)
            A_frequenciesModel = 'fix '

    def pick_GR(e):
        if GR_menu.get() == " ":
            GR_entry_1.delete(first=0, last=1000)
            GR_entry_1.configure(state="disabled")
            GR_entry_2.delete(first=0, last=1000)
            GR_entry_2.configure(state="disabled")
        elif GR_menu.get() == "0":
            GR_entry_1.delete(first=0, last=1000)
            GR_entry_1.configure(state="normal")
            GR_entry_2.delete(first=0, last=1000)
            GR_entry_2.configure(state="disabled")
        elif GR_menu.get() == "2":
            GR_entry_1.delete(first=0, last=1000)
            GR_entry_1.configure(state="disabled")
            GR_entry_2.delete(first=0, last=1000)
            GR_entry_2.configure(state="normal")

    def pick_G(e):
        if G_menu.get() == "fix":
            G_menu_1.delete(first=0,last=1000)
            G_menu_2.current(0)
            G_menu_2.configure(state="disabled")
            G_menu_3.delete(first=0,last=1000)
        if G_menu.get() == "uniform":
            G_menu_1.delete(first=0,last=1000)
            G_menu_2.current(0)
            G_menu_2.configure(state="disabled")
            G_menu_3.delete(first=0,last=1000)
        if G_menu.get() == "norm":
            G_menu_1.delete(first=0,last=1000)
            G_menu_2.current(0)
            G_menu_2.configure(state="enabled")
            G_menu_3.delete(first=0,last=1000)
        if G_menu.get() == "exp":
            G_menu_1.delete(first=0,last=1000)
            G_menu_2.current(0)
            G_menu_2.configure(state="enabled")
            G_menu_3.delete(first=0,last=1000)
        if G_menu.get() == "gamma":
            G_menu_1.delete(first=0,last=1000)
            G_menu_2.current(0)
            G_menu_2.configure(state="enabled")
            G_menu_3.delete(first=0,last=1000)
        if G_menu.get() == "beta":
            G_menu_1.delete(first=0,last=1000)
            G_menu_2.current(0)
            G_menu_2.configure(state="enabled")
            G_menu_3.delete(first=0,last=1000)
        if G_menu.get() == "dirichlet":
            G_menu_1.delete(first=0,last=1000)
            G_menu_2.current(0)
            G_menu_2.configure(state="disabled")
            G_menu_3.delete(first=0,last=1000)
            G_menu_3.configure(state="disabled")
    
    def pick_I(e):
        if I_menu.get() == "fix":
            I_menu_1.delete(first=0,last=1000)
            I_menu_2.current(0)
            I_menu_2.configure(state="disabled")
            I_menu_3.delete(first=0,last=1000)
        if I_menu.get() == "uniform":
            I_menu_1.delete(first=0,last=1000)
            I_menu_2.current(0)
            I_menu_2.configure(state="disabled")
            I_menu_3.delete(first=0,last=1000)
        if I_menu.get() == "norm":
            I_menu_1.delete(first=0,last=1000)
            I_menu_2.current(0)
            I_menu_2.configure(state="enabled")
            I_menu_3.delete(first=0,last=1000)
        if I_menu.get() == "exp":
            I_menu_1.delete(first=0,last=1000)
            I_menu_2.current(0)
            I_menu_2.configure(state="enabled")
            I_menu_3.delete(first=0,last=1000)
        if I_menu.get() == "gamma":
            I_menu_1.delete(first=0,last=1000)
            I_menu_2.current(0)
            I_menu_2.configure(state="enabled")
            I_menu_3.delete(first=0,last=1000)
        if I_menu.get() == "beta":
            I_menu_1.delete(first=0,last=1000)
            I_menu_2.current(0)
            I_menu_2.configure(state="enabled")
            I_menu_3.delete(first=0,last=1000)
        if I_menu.get() == "dirichlet":
            I_menu_1.delete(first=0,last=1000)
            I_menu_2.current(0)
            I_menu_2.configure(state="disabled")
            I_menu_3.delete(first=0,last=1000)
            I_menu_3.configure(state="disabled")

    def Hide_General_Settings():
        global is_on
        # determin is on ot off
        if is_on:
            Additional_parameters.config(text="- Additional parameters")
            is_on = False
            #Oplab1.grid(row=500, column=1, sticky="W", padx=5)
            #Oplab2.grid(row=501, column=1, sticky="W", padx=10)
            #Numb_proce.grid(row=502, column=1, sticky='W', padx=150, pady=5)
            #Intt_Label.grid(row=502, column=1, sticky='W', padx=120)
            Oplab3.grid(row=503, column=1, sticky="W", padx=5)
            Oplab4.grid(row=504, column=1, sticky="W", padx=10)
            sampling.grid(row=505, column=1, sticky='W', padx=150, pady=5)
            Oplab5.grid(row=506, column=1, sticky="W", padx=5)
            G_time.grid(row=507, column=1, sticky='W', padx=150, pady=5)
            Oplab6.grid(row=508, column=1, sticky="W", padx=5)
            Oplab7.grid(row=509, column=1, sticky="W", padx=10)
            G_menu.grid(row=510, column=1, sticky='W', padx=10, pady=5)
            Oplab8.grid(row=509, column=1, sticky="W", padx=10)
            G_menu_1.grid(row=510, column=1, sticky='W', padx=140, pady=5)
            Oplab9.grid(row=509, column=1, sticky="W", padx=285)
            G_menu_2.grid(row=510, column=1, sticky='W', padx=285, pady=5)
            Oplab10.grid(row=509, column=1, sticky="E", padx=140)
            G_menu_3.grid(row=510, column=1, sticky='E', padx=140, pady=5)
            Oplab11.grid(row=517, column=1, sticky="W", padx=5)
            Oplab12.grid(row=518, column=1, sticky="W", padx=10)
            I_menu.grid(row=519, column=1, sticky='W', padx=10, pady=5)
            Oplab13.grid(row=518, column=1, sticky="W", padx=10)
            I_menu_1.grid(row=519, column=1, sticky='W', padx=140)
            Oplab14.grid(row=518, column=1, sticky="W", padx=285, pady=5)
            I_menu_2.grid(row=519, column=1, sticky='W', padx=285)
            Oplab15.grid(row=518, column=1, sticky="E", padx=140)
            I_menu_3.grid(row=519, column=1, sticky='E', padx=140, pady=5)
            Oplab16.grid(row=520, column=1, sticky='W', padx=5)
            Oplab17.grid(row=521, column=1, sticky='W', padx=10)
            GR_menu.grid(row=522, column=1, sticky='W', padx=150, pady=5)
            Oplab18.grid(row=523, column=1, sticky='W', padx=5)
            Oplab19.grid(row=524, column=1, sticky='W', padx=10)
            GR_entry_1.grid(row=525, column=1, sticky='W', padx=30, pady=5)
            Oplab20.grid(row=523, column=1, sticky='E', padx=130)
            Oplab21.grid(row=524, column=1, sticky='E', padx=130)
            GR_entry_2.grid(row=525, column=1, sticky='E', padx=160, pady=5)
            Oplab22.grid(row=526, column=1, sticky="W", padx=5)
            Oplab23.grid(row=527, column=1, sticky="W", padx=10)
            MM_entry.grid(row=528, column=1, sticky="W", padx=150, pady=5)
            Oplab24.grid(row=529, column=1, sticky="W", padx=5)
            Oplab25.grid(row=530, column=1, sticky="W", padx=10)
            MR_entry.grid(row=531, column=1, sticky="W", padx=150, pady=5)
            Oplab26.grid(row=532, column=1, sticky="W", padx=5)
            Oplab27.grid(row=533, column=1, sticky="W", padx=10)
            CD_entry.grid(row=534, column=1, sticky="W", padx=150, pady=5)
            updateScrollRegion()

        else:
            Additional_parameters.config(text="+ Additional parameters")
            is_on = True
            #Oplab1.grid_forget()
            #Oplab2.grid_forget()
            #Numb_proce.grid_forget()
            #ntt_Label.grid_forget()
            Oplab3.grid_forget()
            Oplab4.grid_forget()
            sampling.grid_forget()
            Oplab5.grid_forget()
            G_time.grid_forget()
            Oplab6.grid_forget()
            Oplab7.grid_forget()
            G_menu.grid_forget()
            Oplab8.grid_forget()
            G_menu_1.grid_forget()
            Oplab9.grid_forget()
            G_menu_2.grid_forget()
            Oplab10.grid_forget()
            G_menu_3.grid_forget()
            Oplab11.grid_forget()
            Oplab12.grid_forget()
            I_menu.grid_forget()
            Oplab13.grid_forget()
            I_menu_1.grid_forget()
            Oplab14.grid_forget()
            I_menu_2.grid_forget()
            Oplab15.grid_forget()
            I_menu_3.grid_forget()
            Oplab16.grid_forget()
            Oplab17.grid_forget()
            GR_menu.grid_forget()
            Oplab18.grid_forget()
            Oplab19.grid_forget()
            GR_entry_1.grid_forget()
            Oplab20.grid_forget()
            Oplab21.grid_forget()
            GR_entry_2.grid_forget()
            Oplab22.grid_forget()
            Oplab23.grid_forget()
            MM_entry.grid_forget()
            Oplab24.grid_forget()
            Oplab25.grid_forget()
            MR_entry.grid_forget()
            Oplab26.grid_forget()
            Oplab27.grid_forget()
            CD_entry.grid_forget()
            updateScrollRegion()


    def save_selected_values():
        ## Filename check
        try:
            alignment_directory2
            #sys.tracebacklimit = 0
        except (NameError, ValueError):
            popup("ERROR", "Error in name of input file with MSA\nYou must upload an alignment file")

            #sys.tracebacklimit = 0
            #raise ValueError()
        #except NameError:
            #popup("ERROR", "Error in name of input file with MSA\nYou must upload an alignment file")
            #sys.tracebacklimit = 0
            #raise ValueError()

        ## Number of simulations check
        Numb_simu2 = Numb_simu.get()
        if Numb_simu.get() != '':
            if Numb_simu2.isnumeric() == True:
                if int(Numb_simu2) < 100000:
                    if int(Numb_simu2) > 0:
                        lab4_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                        lab4_2.grid(row=6, column=1, sticky="W", padx=9)
                    else:
                        lab4_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                        lab4_2.grid(row=6, column=1, sticky="W", padx=9)
                        popup("ERROR", "Error in number of simulations\nYou must enter a value between 0-100000")
                        sys.tracebacklimit = 0
                        raise ValueError()
                else:
                    lab4_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab4_2.grid(row=6, column=1, sticky="W", padx=9)
                    popup("ERROR", "Error in number of simulations\nYou must enter a value between 0-100000")
                    sys.tracebacklimit = 0
                    raise ValueError()
            else:
                lab4_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab4_2.grid(row=6, column=1, sticky="W", padx=9)
                popup("ERROR", "Error in number of simulations\nYou must enter a number")
                sys.tracebacklimit = 0
                raise ValueError()
        else:
            lab4_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab4_2.grid(row=6, column=1, sticky="W", padx=9)
            popup("ERROR", "Error in number of simulations\nYou must enter the number of simulations")
            sys.tracebacklimit = 0
            raise ValueError()

        ## Haploid/Diploid check
        if str(haploid_menu.get()) != '':
            lab11_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            lab11_2.grid(row=17, column=1, sticky="W", padx=9)
        else:
            lab11_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab11_2.grid(row=17, column=1, sticky="W", padx=9)
            popup("ERROR", "Error in Haploid/Diploid simulated data\nYou must select the haploid or diploid option")
            sys.tracebacklimit = 0
            raise ValueError()

        ## Efective population size check
        Ne2 = Ne.get()
        if Ne.get() != '':
            if Ne2.isnumeric() == True:
                lab12_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                lab12_2.grid(row=19, column=1, sticky="W", padx=9)
            else:
                lab12_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab12_2.grid(row=19, column=1, sticky="W", padx=9)
                popup("ERROR", "Error in effective population size\nYou must enter a number")
                sys.tracebacklimit = 0
                raise ValueError()
        else:
            lab12_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab12_2.grid(row=19, column=1, sticky="W", padx=9)
            popup("ERROR", "Error in effective population size\nYou must enter the number of simulations")
            sys.tracebacklimit = 0
            raise ValueError()

        ## Aminoacid substitution rate per site
        if str(SR_menu.get()) == '':
            lab20_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab20_2.grid(row=31, column=1, sticky="W", padx=9)
            popup("ERROR", "Error in Aminoacid Substitution Ratel\nYou must select at least one model")
            sys.tracebacklimit = 0
            raise ValueError()
        else:
            if str(SR_menu.get()) == 'fix' and str(SR_menu_1.get()) == '':
                lab20_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab20_2.grid(row=31, column=1, sticky="W", padx=9)
                popup("ERROR", "Error in Aminoacid Substitution Ratel\nYou must enter a value")
                sys.tracebacklimit = 0
                raise ValueError()
            elif str(SR_menu.get()) == 'norm' and str(SR_menu_2.get()) == 't' and str(SR_menu_3.get()) =='':
                lab20_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab20_2.grid(row=31, column=1, sticky="W", padx=9)
                popup("ERROR", "Error in Aminoacid Substitution Ratel\nIf you select truncated option you must enter a value")
                sys.tracebacklimit = 0
                raise ValueError()
            elif str(SR_menu.get()) == 'exp' and str(SR_menu_2.get()) == 't' and str(SR_menu_3.get()) =='':
                lab20_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab20_2.grid(row=31, column=1, sticky="W", padx=9)
                popup("ERROR", "Error in Aminoacid Substitution Ratel\nIf you select truncated option you must enter a value")
                sys.tracebacklimit = 0
                raise ValueError()
            elif str(SR_menu.get()) == 'gamma' and str(SR_menu_2.get()) == 't' and str(SR_menu_3.get()) =='':
                lab20_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab20_2.grid(row=31, column=1, sticky="W", padx=9)
                popup("ERROR", "Error in Aminoacid Substitution Ratel\nIf you select truncated option you must enter a value")
                sys.tracebacklimit = 0
                raise ValueError()
            elif str(SR_menu.get()) == 'beta' and str(SR_menu_2.get()) == 't' and str(SR_menu_3.get()) =='':
                lab20_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab20_2.grid(row=31, column=1, sticky="W", padx=9)
                popup("ERROR", "Error in Aminoacid Substitution Ratel\nIf you select truncated option you must enter a value")
                sys.tracebacklimit = 0
                raise ValueError()
            else:
                lab20_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                lab20_2.grid(row=31, column=1, sticky="W", padx=9)

        ## SubstitutionModel
        if S_model1 != "" or SS_model1 != "":
            lab23_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            lab23_2.grid(row=39, column=1, sticky="W", padx=9)
        else:
            lab23_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab23_2.grid(row=39, column=1, sticky="W", padx=9)
            popup("ERROR", "Error in Substitution Model\nYou must select at least one substitution model")
            sys.tracebacklimit = 0
            raise ValueError()

        ## Amino acid frequencies
        A_frequencies2= A_frequencies.get()
        if A_frequencies2 != "":
            if len(A_frequencies2.split(" ")) != 20:
                lab25_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab25_2.grid(row=43, column=1, sticky="W", padx=9)
                popup("ERROR", "Error in amino acid frequencies\nYou must enter 20 values separated by spaces\nand you have entered " + str(len(A_frequencies2.split(" "))))
                sys.tracebacklimit = 0
                raise ValueError()
            else:
                for e in A_frequencies2.split(" "):
                    if e.isnumeric() == True:
                        if int(e) <=1 :
                            lab25_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                            lab25_2.grid(row=43, column=1, sticky="W", padx=9)
                        else:
                            lab25_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                            lab25_2.grid(row=43, column=1, sticky="W", padx=9)
                            popup("ERROR", "Error in amino acid frequencies\nYou must enter values less or equal to 1")
                            sys.tracebacklimit = 0
                            raise ValueError()
                    else:
                        if e.lstrip("-").isnumeric() == True:
                            lab25_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                            lab25_2.grid(row=43, column=1, sticky="W", padx=9)
                            popup("ERROR", "Error in amino acid frequencies\nYou must enter only positive numbers")
                            sys.tracebacklimit = 0
                            raise ValueError()
                        elif e.replace('.', '', 1).isdigit() == True:
                            lab25_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                            lab25_2.grid(row=43, column=1, sticky="W", padx=9)
                        else:
                            lab25_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                            lab25_2.grid(row=43, column=1, sticky="W", padx=9)
                            popup("ERROR", "Error in amino acid frequencies\nYou must enter only numbers")
                            sys.tracebacklimit = 0
                            raise ValueError()

        else:
            lab25_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab25_2.grid(row=43, column=1, sticky="W", padx=9)
            popup("ERROR", "Error in amino acid frequencies\nYou must enter 20 values separated by spaces")
            sys.tracebacklimit = 0
            raise ValueError()

        ## ABC iteractions
        ABC_itera2 = ABC_itera.get()
        if ABC_itera2 != "":
            if ABC_itera2.isnumeric() == True:
                if int(ABC_itera2) <= int(Numb_simu2):
                    if int(ABC_itera2) > 0:
                        lab27_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                        lab27_2.grid(row=47, column=1, sticky="W", padx=9)
                    else:
                        lab27_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                        lab27_2.grid(row=47, column=1, sticky="W", padx=9)
                        popup("ERROR","Error ABC iterations\nYou must enter a value greater than 0")
                        sys.tracebacklimit = 0
                        raise ValueError()
                else:
                    lab27_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab27_2.grid(row=47, column=1, sticky="W", padx=9)
                    popup("ERROR", "Error ABC iterations\nYou must enter a value less or equal\n to the number of simulations")
                    sys.tracebacklimit = 0
                    raise ValueError()
            else:
                if ABC_itera2.lstrip("-").isnumeric() == True:
                    lab27_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab27_2.grid(row=47, column=1, sticky="W", padx=9)
                    popup("ERROR", "Error ABC iterations\nYou must enter a positive value")
                    sys.tracebacklimit = 0
                    raise ValueError()
                elif ABC_itera2.replace('.', '', 1).isdigit() == True:
                    lab27_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab27_2.grid(row=47, column=1, sticky="W", padx=9)
                    popup("ERROR", "Error ABC iterations\nYou must enter a int number")
                    sys.tracebacklimit = 0
                    raise ValueError()
                else:
                    lab27_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab27_2.grid(row=47, column=1, sticky="W", padx=9)
                    popup("ERROR", "Error ABC iterations\nYou must enter a number")
                    sys.tracebacklimit = 0
                    raise ValueError()

        else:
            lab27_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab27_2.grid(row=47, column=1, sticky="W", padx=9)
            popup("ERROR", "Error ABC iterations\nYou must enter a value")
            sys.tracebacklimit = 0
            raise ValueError()

        ## ABC tolerance
        ABC_tolerance2 = ABC_tolerance.get()
        if ABC_tolerance2 != "":
            if ABC_tolerance2.isnumeric() == True:
                if float(ABC_tolerance2) <= 1:
                    if float(ABC_tolerance2) > 0:
                        lab29_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                        lab29_2.grid(row=50, column=1, sticky="W", padx=9)
                    else:
                        lab29_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                        lab29_2.grid(row=50, column=1, sticky="W", padx=9)
                        popup("ERROR","Error ABC tolerance\nYou must enter a value greater than 0")
                        sys.tracebacklimit = 0
                        raise ValueError()
                else:
                    lab29_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab29_2.grid(row=50, column=1, sticky="W", padx=9)
                    popup("ERROR", "Error ABC tolerance\nYou must enter a value less than 1")
                    sys.tracebacklimit = 0
                    raise ValueError()
            else:
                if ABC_tolerance2.lstrip("-").isnumeric() == True:
                    lab29_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab29_2.grid(row=50, column=1, sticky="W", padx=9)
                    popup("ERROR", "Error ABC tolerance\nYou must enter a positive value")
                    sys.tracebacklimit = 0
                    raise ValueError()
                elif ABC_tolerance2.replace('.', '', 1).isdigit() == True:
                    lab29_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                    lab29_2.grid(row=50, column=1, sticky="W", padx=9)
                else:
                    lab29_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab29_2.grid(row=50, column=1, sticky="W", padx=9)
                    popup("ERROR", "Error ABC tolerance\nYou must enter a number")
                    sys.tracebacklimit = 0
                    raise ValueError()

        else:
            lab29_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab29_2.grid(row=50, column=1, sticky="W", padx=9)
            popup("ERROR", "Error ABC tolerance\nYou must enter a value")
            sys.tracebacklimit = 0
            raise ValueError()

        #ABC Method
        if str(ABCMethod_menu.get()) != '':
            lab31_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            lab31_2.grid(row=53, column=1, sticky="W", padx=9)
        else:
            lab31_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab31_2.grid(row=53, column=1, sticky="W", padx=9)
            popup("ERROR", "Error in ABC method\nYou must select ABC method option")
            sys.tracebacklimit = 0
            raise ValueError()

        ## Summary stadistics
        if SStadistics != "":
            lab39_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            lab39_2.grid(row=64, column=1, sticky="W", padx=9)
        else:
            lab39_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab39_2.grid(row=64, column=1, sticky="W", padx=9)
            popup("ERROR", "Error in Summary Stadistics\nYou must select at least 1 summary stadistics")
            sys.tracebacklimit = 0
            raise ValueError()

        ## Multipages
        if str(MultiPages_menu.get()) != "":
            lab42_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            lab42_2.grid(row=64, column=1, sticky="W", padx=9)
        else:
            lab42_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab42_2.grid(row=64, column=1, sticky="W", padx=9)
            popup("ERROR", "Error in multiple pages\nYou must select a multiple pages value")
            sys.tracebacklimit = 0
            raise ValueError()


        file = open('Settings.txt', 'w')
        file.write('##########################################################################################################################' + '\n')
        file.write('##### Settings file for ProtModel' + '\n')
        file.write('##### Selection of the best-fitting empirical or structural substitution model' + '\n')
        file.write('##### David Ferreiro ' + '\n')
        file.write('##### (c) 2022' + '\n')
        file.write('##### Contact: david.ferreiro.garcia@uvigo.es / ferreirogarciadavid@gmail.com' + '\n')
        file.write('##### ' + '\n')
        file.write('##### Parameters with an "*" are mandatory (need to be specified)' + '\n')
        file.write('##### Parameter values must be introduced immediately after the "="' '\n')
        file.write('##########################################################################################################################' + '\n' + '\n')
        file.write('#########################################' + '\n')
        file.write('### Settings for the simulation phase ###' + '\n')
        file.write('#########################################' + '\n')
        file.write('### Prior distributions available: fix #, uniform # #, gamma # # (t) # #, beta # # (t) # #, normal # # (t) # #, exponential # (t) # #, dirichlet n#' + '\n')
        file.write('### (t) indicates that the distribution can be truncated through the following two (lowest highest) values.' + '\n' + '\n' + '\n')
        file.write('### Target alignment file ###		# phylip format, see documentation for details' + '\n')
        file.write('*NameOfPhylipFile=' + str(alignment_directory2) + '\n' + '\n' + '\n')
        file.write('### Total number of simulations ###' + '\n')
        file.write('*NumberOfSimulations=' + Numb_simu.get() + '\n' + '\n')
        file.write('# Consideration of indels. Ignored, or are considered as a New State' + '\n')
        file.write('*Indels=' + indells_menu.get() + '\n' + '\n')
        file.write('# Number of available processors to run the simulations in parallel. 1 by default' + '\n')
        file.write('*NumberOfProcessors=' + Numb_proce.get() + '\n' + '\n')
        file.write('# Save simulated data. No, Yes (requires space in the disk)' + '\n')
        file.write('*SaveSimulations=' + Save_S_menu.get() + '\n' + '\n')
        file.write('# Show running information on the screen. No, Yes (it increase the running time)' + '\n')
        file.write('*ShowInformationScreen=' + Show_RI_menu.get() + '\n' + '\n' )
        file.write('### Demographic settings ###' + '\n' + '\n')
        file.write('# Haploid or Diploid data (haploid=1, diploid=2)' + '\n')

        if haploid_menu.get() == 'Haploid':
            haploid_value = 1
        elif haploid_menu.get() == 'Diploid':
            haploid_value = 2

        file.write('*Haploid/Diploid=' + str(haploid_value) + '\n' + '\n')
        file.write('# Population size (i.e., 1000)' + '\n')
        file.write('*PopulationSize=' + Ne.get() + '\n' + '\n')
        file.write('# Logitudinal sampling. Requires GenerationTime. See documentation for details' + '\n')
        file.write('DatedTips=' + sampling.get() + '\n' + '\n')
        file.write('# Generation time. fix, uniform; i.e., uniform 500 1000' + '\n')
        file.write('GenerationTime=' + G_time.get() + '\n')
        file.write('# Exponential growth rate or Demographic periods. See documentation for details'+ '\n')
        file.write('GrowthRate=' + str(GR_menu.get()) + ' ' + str(GR_entry_1.get()) + str(GR_entry_2.get()) + '\n')
        file.write('# Migration model and population structure. See documentation for details'+ '\n')
        file.write('MigrationModel=' + MM_entry.get() + '\n')
        file.write('# Migration rate (constant or variable with time according to temporal periods). See documentation for details'+ '\n')
        file.write('MigrationRate=' + MR_entry.get() + '\n')
        file.write('# Events of convergence of demes. See documentation for details'+ '\n')
        file.write('ConvergenceDemes=' + CD_entry.get() + '\n' + '\n')
        #file.write('# Recombination rate per site. fix, uniform, gamma, beta, normal, exponential; i.e., gamma 0.02 0.5 t 1.3e-07 9.3e-07. -PARAMETER TO BE ESTIMATED-' + '\n')
        #file.write('*RecombinationRate=' + recombination_menu.get() + ' ' + recombination_menu_1.get() + ' ' + recombination_menu_2.get() + recombination_menu_3.get() + '\n' + '\n' + '\n')
        file.write('### Protein substitution model settings ###' + '\n' + '\n')
        file.write('# Amino acid substitution rate. i.e., fix 7.0e-6.' + '\n')
        file.write('*SubstitutionRate=' + SR_menu.get() + ' ' + SR_menu_1.get() + ' ' + SR_menu_2.get() + SR_menu_3.get() + '\n' + '\n')
        file.write('# Model of amino acid substitution (i.e., Blosum62, CpRev, Dayhoff, DayhoffDCMUT, HIVb, HIVw, JTT, JonesDCMUT, LG, Mtart, Mtmam, Mtrev24, RtRev, VT, WAG, UserEAAM)' + '\n')
        file.write('*SubstitutionModel=' + str(S_model1) + '\n' + '\n')
        file.write('# Structural amino acid substitution (Fitness and Neutral)' + '\n')
        file.write('*StructuralSubstitutionModel=' + str(SS_model1) + '\n' + '\n')
        file.write('# Amino acid frequencies. fix or dirichlet. By default equally distributed frequencies. i.e., dirichlet 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1' + '\n')
        file.write('*AminoacidFrequencies='+ str(A_frequenciesModel) + A_frequencies.get() + '\n' + '\n')
        file.write('# Rate of heteregeneity across sites, +G. fix, uniform, gamma, beta, normal, exponential; i.e., fix 0.6' + '\n')
        file.write('RateHetSites=' + G_menu.get() + ' ' + G_menu_1.get() + ' ' + G_menu_2.get() + G_menu_3.get() + '\n' + '\n')
        file.write('# Proportion of invariable sites, +I. fix, uniform, gamma, beta, normal, exponential; i.e., exponential 0.002 t 0 1.0' + '\n')
        file.write('PropInvSites=' + I_menu.get() + ' ' + I_menu_1.get() + ' ' + I_menu_2.get() + I_menu_3.get() + '\n' )
        file.write('# PDB protein structure used to structural substitution models. See documentation for details'+ '\n')
        file.write('*Template=' + str(pdb_file1) + '\n')
        file.write('# Chain of the PDB protein structure used to structural substitution models. See documentation for details' + '\n')
        file.write('*Chain=' + str(chain1) + '\n')
        file.write('# GMRCA input file. See documentation for details'+ '\n')
        file.write('GMRCA=' + str(GMRCA1) + '\n')
        file.write('#########################################' + '\n')
        file.write('### Settings for the estimation phase ###' + '\n')
        file.write('#########################################' + '\n' + '\n')
        file.write('### ABC settings ###' + '\n')
        file.write('#ABC iterations. Number of simulations to consider (Iterations <= NumberOfSimulations)' + '\n')
        file.write('*ABCIterations=' + ABC_itera2 + '\n' + '\n')
        file.write('#ABC tolerance. Number of simulations clostest to real data to retain in the ABC procedure (Tolerance < NumberOfSimulations)' + '\n')
        file.write('*ABCTolerance=' + ABC_tolerance2 + '\n' + '\n')
        file.write('#ABC method (rejection, loclinear). See documentation for details' + '\n')
        file.write('*ABCMethod=' + str(ABCMethod_menu.get()) + '\n' + '\n')
        file.write('#Summary statistics to use. See documentation for details' + '\n')
        file.write('*SummaryStatistics=' + SStadistics + '\n' + '\n')
        file.write('### Graphical settings ###' + '\n'  + '\n')
        file.write('#Multiple pages. PDF documents with multiple pages (no=0, yes=1)' + '\n')
        file.write('*MultiPage=' + str(MultiPages_menu.get()) + '\n' + '\n')
        file.close()

        # Toplevel object which will be treated as a new window
        Running_window = Toplevel(second_frame)
        # sets the title of the Toplevel widget
        Running_window.title("Running ProtModel")
        Running_window.config(bg="gray90")
        # sets the geometry of toplevel
        Running_window.geometry("480x480")
        Running_window.geometry(f"480x480+{root.winfo_x()+110}+{root.winfo_y()-150}")
        # Running_window.resizable(False, False)
        # create a frame
        Runningframe = ttk.Frame(Running_window, padding=(5, 150, 5, 5))
        Runningframe.grid(column=0, row=0, sticky=(N, S, E, W))

        lab_r_0 = Label(Running_window, text="", font="Arial 10 bold", bg="gray90")
        lab_r_0.grid(row=0, column=1, columnspan=2)
        lab_r_1 = Label(Running_window, text="Progress", font=('Avenir Next', 14, 'bold'), bg="gray90")
        lab_r_1.grid(row=1, column=1, sticky=E, pady=5)
        s = ttk.Style()
        s.theme_use('clam')
        s.configure("red.Horizontal.TProgressbar", troughcolor='white',
                    bordercolor='darkblue', background='lightblue', lightcolor='darkblue',
                    darkcolor='darkblue')
        Progress = ttk.Progressbar(Running_window, style='red.Horizontal.TProgressbar', orient=HORIZONTAL,
                                   length=426, mode='indeterminate')
        Progress.grid(row=2, column=1, columnspan=2, padx=27)
        lab_r_2 = Label(Running_window, text="", font=('Avenir Next', 14, 'bold'), bg="gray90")
        lab_r_2.grid(row=1, column=2, sticky=W)
        #Progress.start(20)
        #Progress['value'] = 0
        lab_r_2.configure(text="0%")
        lab_r_3 = Label(Running_window, text="", font=('Avenir Next', 12), bg="gray90")
        lab_r_3.place(x=24, y=75)

        text = "S"

        start_button = Button()

        loading = ' ,.,..,...'.split(',')

        def call():
            global w  # index number
            if Progress['value'] != 100:  # if it is in the range of the list
                lab_r_3.config(text="Loading " + loading[w],
                               font=('Avenir Next', 12))  # change the text to the corresponding index element
                w += 1  # increase index number
                Running_window.after(500, call)  # repeat the function after 2 seconds
                if w == 4:
                    w = 0
            else:
                lab_r_3.config(text="Finish",
                               font=('Avenir Next', 12))  # change the text to the corresponding index element
                #print('Done')  # if out of index range, then dont repeat
                #popup("DONE!!\n, ProtModel has finished correctly ")
                
        def reanalyze():
            pass

        #os.system('chmod +x ProtModelGeneral.py')
        text = Text(Running_window, height=25, width=60, highlightthickness=1, highlightbackground='black',borderwidth=10)
        text.grid(row=3, column=1, columnspan=2, padx=15, pady=35)
        
        def actualizartext():
            threading.Timer(1.0, actualizartext).start()
            
        def start_ProtModel():
            Progress['value'] = 0
            lab_r_2.configure(text="0%")
            call()
            Progress.start(20)
            #threading.Timer(1.0, start_ProtModel).start()
            #global p
            #p = sub.Popen('./ProtModel-1.py', stdout=sub.PIPE, stderr=sub.PIPE, bufsize=0) #sin stderr
            p = sub.Popen('./ProtModelGeneral-M.py', stdout=sub.PIPE, bufsize=0, universal_newlines=True) #bufsize=1 sin stderr
            #output = p.communicate()
            #stdout = output[0]
            #stderr = output[1]
            
            #text.insert(END, stdout)
            #text.insert(END, stderr)
            #text.yview_pickplace("end")
            #actualizartext()
            def stopsubprocess():
                p.terminate()
                Progress.stop()
                popup("STOPED\n", "ProtModel has been stopped")
                lab_r_2.configure(text="0%")
                Progress.config(mode='determinate')
                Progress['value']=0
                stopbutton.destroy()
                

            stopbutton = Button(Running_window, text="Stop", font=('Avenir Next', 12, 'bold'), command=stopsubprocess,
                                width=50, relief='flat', focuscolor='red', bg='white', fg='black',
                                borderless=1)  # command=start, relief='sunken', bg='white', bordercolor='darkblue', borderless=1, focuscolor='', width=50)
            stopbutton.place(x=390, y=80)


            while True:
                realtime_output = p.stdout.readline()

                #realtime_error = p.stderr.readline()



                #if realtime_output.decode("utf-8") != "Finish!":
                    #print(realtime_output.decode("utf-8"))
                #if realtime_output == 'Starting simulations based on empirical amino acid substitution models':
                    #Progress['value'] = 20
                    #lab_r_2.configure(text="20%")

                if realtime_output:
                    #text = Text(Running_window, height=25, width=60, highlightthickness=1, highlightbackground='black',borderwidth=10)
                    #text.grid(row=3, column=1, columnspan=2, padx=15, pady=35)
                    text.insert("end", realtime_output)
                    text.yview_pickplace("end")
                    t = realtime_output#.decode("utf-8")
                    if t.startswith('Alignment file exists'):
                        lab_r_2.configure(text="5%")
                    if t.startswith('Starting simulations based on ' + str(S_model1) + ' amino acid substitution model'):
                        lab_r_2.configure(text="10%")
                    if t.startswith('Starting simulations based on ' + str(SS_model[0]) + ' amino acid substitution model'):
                        lab_r_2.configure(text="25%")
                    if t.startswith('Starting simulations based on ' + str(SS_model[1]) + ' amino acid substitution model'):
                        lab_r_2.configure(text="50%")
                    if t.startswith('Summary statistics calculation'):
                        lab_r_2.configure(text="75%")

                else:
                    if Progress['value'] != 0:
                        Progress.stop()
                        popup2("DONE!!\n", "ProtModel has finished correctly")
                        lab_r_2.configure(text="100%")
                        Progress.config(mode='determinate')
                        Progress['value']=100
                        break
                    
            if Progress['value']==100:
                #reanalyse=Button(Running_window, text="Re-analyse", font=('Avenir Next', 12, 'bold'), width=80,relief='flat', focuscolor='darkblue', bg='lightblue', fg='black',borderless=1)  # command=start, relief='sunken', bg='white', bordercolor='darkblue', borderless=1, focuscolor='', width=50)
                #reanalyse.place(x=290, y=80)
                stopbutton.destroy()


                




            #output = p.communicate()
            #stdout = output[0]
            #stderr = output[1]

            #if stderr.decode("utf-8") == "":
                #pass
            #else:
                #stderr = '\n\n*************************************************************\n' + '*************************************************************\n' + '*************************************************************\n\n' + stderr.decode("utf-8") + '\n*************************************************************\n' + '*************************************************************\n' + '*************************************************************'

            #text = Text(Running_window, height=25, width=60, highlightthickness=1, highlightbackground='black', borderwidth=10)
            #text.grid(row=3, column=1, columnspan=2, padx=15, pady=35)
            #text.insert(END, stdout)
            #text.insert(END, stderr)
            #text.yview_pickplace("end")

            #value = text.get("end-62c", "end-1c")
            #if value.find("*************************************************************") == -1: #-1 means not found
                #Progress['value'] = 0
                #lab_r_2.configure(text="0%")
            #else:
                #Progress['value'] = 100
                #lab_r_2.configure(text="100%")
                #call()  # call the function initially

        # p = sub.Popen('./ProtModel-2.py',stdout=sub.PIPE,stderr=sub.PIPE)
        # output = p.communicate()
        # stdout=output[0]
        # stderr=output[1]
        # if stderr.decode("utf-8") == "":
        #    pass
        # else:
        #    stderr='*************************************************************\n' + stderr.decode("utf-8") + '*************************************************************\n'

        # text.grid(row=3, column=1, columnspan=2, padx=15, pady=35)
        # text.insert(END, stdout)
        # text.insert(END, stderr)
        # text.yview_pickplace("end")

        # Progress['value'] = 50
        # lab_r_2.configure(text="50%")

        # p = sub.Popen('./ProtModel-3.py',stdout=sub.PIPE,stderr=sub.PIPE)
        # output = p.communicate()
        # stdout=output[0]
        # stderr=output[1]
        # if stderr.decode("utf-8") == "":
        #    pass
        # else:
        #    stderr='*************************************************************\n' + stderr.decode("utf-8") + '*************************************************************\n'

        # text.grid(row=3, column=1, columnspan=2, padx=15, pady=35)
        # text.insert(END, stdout)
        # text.insert(END, stderr)
        # text.yview_pickplace("end")

        # Progress['value'] = 80
        # lab_r_2.configure(text="80%")

        # p = sub.Popen('./ProtModel-4.py',stdout=sub.PIPE,stderr=sub.PIPE)
        # output = p.communicate()
        # stdout=output[0]
        # stderr=output[1]
        # if stderr.decode("utf-8") == "":
        #    pass
        # else:
        #    stderr='*************************************************************\n' + stderr.decode("utf-8") + '*************************************************************\n'

        # text.grid(row=3, column=1, columnspan=2, padx=15, pady=35)
        # text.insert(END, stdout)
        # text.insert(END, stderr)
        # text.yview_pickplace("end")

        # Progress['value'] = 100
        # lab_r_2.configure(text="100%")

        start = Button(Running_window, text="Start", font=('Avenir Next', 12, 'bold'), command=lambda: threading.Thread(target=start_ProtModel).start(), width=50,relief='flat', focuscolor='darkblue', bg='lightblue', fg='black',borderless=1)  # command=start, relief='sunken', bg='white', bordercolor='darkblue', borderless=1, focuscolor='', width=50)
        start.place(x=390, y=80)
        #start.focus()

        Running_window.mainloop()


    #Pmw.initialise(root)
    lab0 = Label(second_frame, text="  -------------------------------------------------------------------------------------------------------------------------------------------",font="Arial 10", anchor='w', fg='gray70', bg='gray90')
    lab0.grid(row=0, column=1, sticky="W")
    lab0.bind('<MouseWheel>', on_vertical)

    lab2 = Label(second_frame, text="   Necessary settings", font=('Arial', 16, 'bold'),bg='gray90')
    lab2.grid(row=1, column=1, sticky="W", padx=5, pady=5)
    lab2.bind('<MouseWheel>', on_vertical)

    #subsection1 = Label(second_frame, text="   ### Settings for the simulation phase ###", font=("Arial", "14", "italic"), anchor='w', bg='gray90')
    #subsection1.grid(row=2, column=1, sticky="W", padx=5, pady=5)
    #subsection1.bind('<MouseWheel>', on_vertical)

    lab3 = Label(second_frame, text="* Name of input file with MSA:", font="Arial 14", bg='gray90')
    lab3.grid(row=3, column=1, sticky="W", padx=5, pady=2)
    lab3.bind('<MouseWheel>', on_vertical)
    lab3_2 = Label(second_frame, text="", font="Arial 14", anchor='w', bg='gray90')
    lab3_2.grid(row=4, column=1, sticky="W", padx=35)
    lab3_2.bind('<MouseWheel>', on_vertical)
    upload = Button(second_frame, text="Upload ↑", command=Name_file, bg='white', borderless=1, focuscolor='')
    upload.grid(row=4, column=1, sticky="E", padx=140)
    upload.bind('<MouseWheel>', on_vertical)

    lab4 = Label(second_frame, text="* Number of simulations (integrer):", font="Arial 14", anchor='w', bg='gray90')
    lab4.grid(row=5, column=1, sticky="W", padx=5)
    lab4.bind('<MouseWheel>', on_vertical)
    lab4_2 = Label(second_frame, text="", font="Arial 14", anchor='w', bg='gray90')
    lab4_2.bind('<MouseWheel>', on_vertical)
    Numb_simu = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    Numb_simu.grid(row=6, column=1, sticky='W', padx=150, pady=5)
    Numb_simu.bind('<MouseWheel>', on_vertical)

    lab5 = Label(second_frame, text="* How to consider indels (gaps):", font="Arial 14", anchor='w', bg='gray90')
    lab5.grid(row=7, column=1, sticky="W", padx=5)
    lab5.bind('<MouseWheel>', on_vertical)
    lab6 = Label(second_frame,text="Ignored(by default), New state (indels are considered as a new state):",font="Arial 10", anchor='w', bg='gray90')
    lab6.grid(row=8, column=1, sticky="W", padx=10)
    lab6.bind('<MouseWheel>', on_vertical)
    indells_menu = ttk.Combobox(second_frame, value=gaps, state="readonly")
    indells_menu.current(0)
    indells_menu.grid(row=9, column=1, sticky='W', padx=150,pady=5)
    indells_menu.bind('<MouseWheel>', on_vertical)

    lab7 = Label(second_frame, text="* Save simulations?:", font="Arial 14", anchor='w', bg='gray90')
    lab7.grid(row=10, column=1, sticky="W", padx=5)
    lab7.bind('<MouseWheel>', on_vertical)
    lab8 = Label(second_frame, text="No(by default), Yes (this option requires space in the disk)", font="Arial 10",anchor='w', bg='gray90')
    lab8.grid(row=11, column=1, sticky="W", padx=10)
    lab8.bind('<MouseWheel>', on_vertical)
    Save_S_menu = ttk.Combobox(second_frame, value=S_Simulations, state="readonly")
    Save_S_menu.current(0)
    Save_S_menu.grid(row=12, column=1, sticky='W', padx=150,pady=5)
    Save_S_menu.bind('<MouseWheel>', on_vertical)

    lab9 = Label(second_frame, text="* Show running information on the screen?:", font="Arial 14", anchor='w',bg='gray90')
    lab9.grid(row=13, column=1, sticky="W", padx=5)
    lab9.bind('<MouseWheel>', on_vertical)
    lab10 = Label(second_frame, text="No(by default), Yes (this option increase the running time)", font="Arial 10",anchor='w', bg='gray90')
    lab10.grid(row=14, column=1, sticky="W", padx=10)
    lab10.bind('<MouseWheel>', on_vertical)
    Show_RI_menu = ttk.Combobox(second_frame, value=S_Running, state="readonly")
    Show_RI_menu.current(0)
    Show_RI_menu.grid(row=15, column=1, sticky='W', padx=150,pady=5)
    Show_RI_menu.bind('<MouseWheel>', on_vertical)
    
    lab11 = Label(second_frame, text="* Haploid or Diploid simulated data", font=('Arial', 14), anchor='w', bg='gray90')
    lab11.grid(row=16, column=1, sticky="W", padx=5)
    lab11.bind('<MouseWheel>', on_vertical)
    lab11_2 = Label(second_frame, text="", font="Arial 14", anchor='w', bg='gray90')
    lab11_2.grid(row=17, column=1, sticky='W', padx=150)
    lab11_2.bind('<MouseWheel>', on_vertical)
    haploid_menu = ttk.Combobox(second_frame, value=haploid, state="readonly")
    haploid_menu.grid(row=17, column=1, sticky='W', padx=150,pady=5)
    haploid_menu.bind('<MouseWheel>', on_vertical)
    
    lab12 = Label(second_frame, text="* Effective population size (N)(integrer)", font=('Arial', 14), anchor='w', bg='gray90')
    lab12.grid(row=18, column=1, sticky="W", padx=5)
    lab12.bind('<MouseWheel>', on_vertical)
    lab12_2 = Label(second_frame, text="", font="Arial 14", anchor='w', bg='gray90')
    lab12_2.grid(row=19, column=1, sticky='W', padx=150)
    lab12_2.bind('<MouseWheel>', on_vertical)
    Ne = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    Ne.grid(row=19, column=1, sticky='W', padx=150,pady=5)
    Ne.bind('<MouseWheel>', on_vertical)

    Oplab1 = Label(second_frame, text="* Number of processors", font="Arial 14", anchor='w', bg='gray90')
    Oplab1.grid(row=20, column=1, sticky="W", padx=5)
    Oplab1.bind('<MouseWheel>', on_vertical)
    Oplab2 = Label(second_frame, text="Max number of proccesors selected by default", font="Arial 10", anchor='w', bg='gray90')
    Oplab2.grid(row=21, column=1, sticky="W", padx=5)
    Oplab2.bind('<MouseWheel>', on_vertical)
    Numb_proce = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    Numb_proce.insert(END, os.cpu_count())
    Numb_proce.grid(row=22, column=1, sticky="W", padx=150, pady=5)
    Numb_proce.bind('<MouseWheel>', on_vertical)


    #lab13 = Label(second_frame, text="* Recombination Rate per site", font=('Arial', 12), anchor='w', bg='gray90')
    #lab13.grid(row=20, column=1, sticky="W", padx=5)
    #lab14 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    #lab14.grid(row=21, column=1, sticky="W", padx=10)
    #recombination_menu = ttk.Combobox(second_frame, value=distribution, state="readonly")
    #recombination_menu.bind("<<ComboboxSelected>>", pick_reco)
    #recombination_menu.grid(row=22, column=1)
    #lab15 = Label(second_frame, text="Example: 1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    #lab15.grid(row=23, column=1, sticky="W", padx=10)
    #recombination_menu_1 = Entry(second_frame,width=20,bd=3)
    #recombination_menu_1.grid(row=24, column=1)
    #lab16 = Label(second_frame, text="Not truncated or truncated (t)", font=('Arial', 10), anchor='w', bg='gray90')
    #lab16.grid(row=25, column=1, sticky="W", padx=10)
    #recombination_menu_2 = ttk.Combobox(second_frame, value=truncated, state="readonly")
    #recombination_menu_2.grid(row=26, column=1)
    #lab17 = Label(second_frame, text="Example: 1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    #lab17.grid(row=27, column=1, sticky="W", padx=10)
    #recombination_menu_3 = Entry(second_frame,width=20,bd=3)
    #recombination_menu_3.grid(row=28, column=1)
    
    lab18 = Label(second_frame, text="* Amino acid substitution rate per site", font=('Arial', 14), anchor='w', bg='gray90')
    lab18.grid(row=29, column=1, sticky="W", padx=5)
    lab18.bind('<MouseWheel>', on_vertical)
    lab19 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab19.grid(row=30, column=1, sticky="W", padx=10)
    lab19.bind('<MouseWheel>', on_vertical)
    SR_menu = ttk.Combobox(second_frame, value=distribution, state="readonly")
    SR_menu.config(width=10)
    SR_menu.bind("<<ComboboxSelected>>", pick_SR)
    SR_menu.grid(row=32, column=1, sticky="W", padx=10,pady=5)
    SR_menu.bind('<MouseWheel>', on_vertical)
    lab20= Label(second_frame, text="          Example: norm              1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    lab20.grid(row=31, column=1, sticky="W", padx=10)
    lab20.bind('<MouseWheel>', on_vertical)
    lab20_2= Label(second_frame, text="", font=('Arial', 10), anchor='w', bg='gray90')
    lab20_2.grid(row=31, column=1, sticky="W", padx=10)
    lab20_2.bind('<MouseWheel>', on_vertical)
    SR_menu_1 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    SR_menu_1.config(width=12)
    SR_menu_1.grid(row=32, column=1, sticky="W", padx=140, pady=5)
    SR_menu_1.bind('<MouseWheel>', on_vertical)
    lab21 = Label(second_frame, text="   (t)", font=('Arial', 10), anchor='w', bg='gray90')
    lab21.grid(row=31, column=1, sticky="W", padx=285)
    lab21.bind('<MouseWheel>', on_vertical)
    SR_menu_2 = ttk.Combobox(second_frame, value=truncated, state="readonly")
    SR_menu_2.config(width=3)
    SR_menu_2.grid(row=32, column=1, sticky="W", padx=285, pady=5)
    SR_menu_2.bind('<MouseWheel>', on_vertical)
    lab22 = Label(second_frame, text="1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    lab22.grid(row=31, column=1, sticky="E", padx=140)
    lab22.bind('<MouseWheel>', on_vertical)
    SR_menu_3 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    SR_menu_3.config(width=12)
    SR_menu_3.grid(row=32, column=1, sticky="E", padx=140, pady=5)
    SR_menu_3.bind('<MouseWheel>', on_vertical)
    
    lab23 = Label(second_frame, text="* Model of amino acid substitution and template selection", font=('Arial', 14), anchor='w', bg='gray90')
    lab23.grid(row=38, column=1, sticky="W", padx=5)
    lab23.bind('<MouseWheel>', on_vertical)
    #S_model_menu = ttk.Combobox(second_frame, value=S_model)
    #S_model_menu.grid(row=39, column=1)
    SubModel = Button(second_frame, text="Substitution Model", command=openSM, bg='white', bordercolor='gray90', borderless=1, focuscolor='', width=150)
    SubModel.grid(row=39, column=1, sticky='W', padx=170, pady=5)
    SubModel.bind('<MouseWheel>', on_vertical)
    lab23_2 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab23_2.grid(row=39, column=1, sticky="W", padx=10, pady=5)
    lab23_2.bind('<MouseWheel>', on_vertical)
    
    lab24 = Label(second_frame, text="* Amino acid frequencies", font=('Arial', 14), anchor='w', bg='gray90')
    lab24.grid(row=40, column=1, sticky="W", padx=5)
    lab24.bind('<MouseWheel>', on_vertical)
    lab25 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab25.grid(row=41, column=1, sticky="W", padx=10)
    lab25.bind('<MouseWheel>', on_vertical)
    lab25_2 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab25_2.grid(row=43, column=1, sticky="W", padx=10)
    lab25_2.bind('<MouseWheel>', on_vertical)
    A_frequencies_menu=ttk.Combobox(second_frame, value=amino_distribution, state="readonly")
    A_frequencies_menu.bind("<<ComboboxSelected>>", pick_amino)
    A_frequencies_menu.current(0)
    A_frequencies_menu.grid(row=42, column=1, sticky='W', padx=150, pady=5)
    A_frequencies_menu.bind('<MouseWheel>', on_vertical)
    A_frequencies = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    A_frequencies.grid(row=43, column=1, sticky='W', padx=150, pady=5)
    A_frequencies.insert(0, "0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05")
    A_frequencies.configure(state="disabled")
    A_frequencies.bind('<MouseWheel>', on_vertical)



    subsection2 = Label(second_frame, text="    Settings for the estimation phase ", font=("Arial", "16", "italic"), anchor='w', bg='gray90')
    subsection2.grid(row=44, column=1, sticky="W", padx=5, pady=5)
    subsection2.bind('<MouseWheel>', on_vertical)

    lab26 = Label(second_frame, text="* ABC iterations", font=('Arial', 14), anchor='w', bg='gray90')
    lab26.grid(row=45, column=1, sticky="W", padx=5)
    lab26.bind('<MouseWheel>', on_vertical)
    lab27 = Label(second_frame, text="Iterations <= NumberOfSimulations", font="Arial 10", anchor='w', bg='gray90')
    lab27.grid(row=46, column=1, sticky="W", padx=10)
    lab27.bind('<MouseWheel>', on_vertical)
    lab27_2 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab27_2.grid(row=47, column=1, sticky="W", padx=10)
    lab27_2.bind('<MouseWheel>', on_vertical)
    ABC_itera = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    ABC_itera.grid(row=47, column=1, sticky='W', padx=150, pady=5)
    ABC_itera.bind('<MouseWheel>', on_vertical)

    lab28 = Label(second_frame, text="* ABC tolerance", font=('Arial', 14), anchor='w', bg='gray90')
    lab28.grid(row=48, column=1, sticky="W", padx=5)
    lab28.bind('<MouseWheel>', on_vertical)
    lab29 = Label(second_frame, text="% of acepted simulations. Example: 0.01", font="Arial 10", anchor='w', bg='gray90')
    lab29.grid(row=49, column=1, sticky="W", padx=10)
    lab29.bind('<MouseWheel>', on_vertical)
    lab29_2 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab29_2.grid(row=50, column=1, sticky="W", padx=10)
    lab29_2.bind('<MouseWheel>', on_vertical)
    ABC_tolerance = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    ABC_tolerance.grid(row=50, column=1, sticky='W', padx=150, pady=5)
    ABC_tolerance.bind('<MouseWheel>', on_vertical)

    lab30 = Label(second_frame, text="* ABC method", font=('Arial', 14), anchor='w', bg='gray90')
    lab30.grid(row=51, column=1, sticky="W", padx=5)
    lab30.bind('<MouseWheel>', on_vertical)
    lab31 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab31.grid(row=52, column=1, sticky="W", padx=10)
    lab31.bind('<MouseWheel>', on_vertical)
    lab31_2 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab31_2.grid(row=53, column=1, sticky="W", padx=10)
    lab31_2.bind('<MouseWheel>', on_vertical)
    ABCMethod_menu = ttk.Combobox(second_frame, value=ABCMethod, state="readonly")
    ABCMethod_menu.grid(row=53, column=1, sticky='W', padx=150, pady=5)
    ABCMethod_menu.bind('<MouseWheel>', on_vertical)

    lab37 = Label(second_frame,text="  -------------------------------------------------------------------------------------------------------------------------------------------",font="Arial 10", anchor='w', fg='gray70', bg='gray90')
    lab37.grid(row=61, column=1, sticky="W", pady=5)
    lab37.bind('<MouseWheel>', on_vertical)

    lab38 = Label(second_frame, text="* Summary statistics to use.", font=('Arial', 14), anchor='w', bg='gray90')
    lab38.grid(row=62, column=1, sticky="W", padx=5)
    lab38.bind('<MouseWheel>', on_vertical)
    lab39 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab39.grid(row=63, column=1, sticky="W", padx=10)
    lab39.bind('<MouseWheel>', on_vertical)
    lab39_2 =(Label(second_frame, text='', font="Arial 10", anchor='w', bg='gray90'))
    lab39_2.grid(row=64, column=1, sticky="W", padx=10, pady=5)
    lab39_2.bind('<MouseWheel>', on_vertical)
    #Summary_menu = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    #Summary_menu.grid(row=64, column=1, sticky='W', padx=150)
    #Summary_menu.bind('<MouseWheel>', on_vertical)
    SubModels = Button(second_frame, text="Sumary stadistics", command=openSS, bg='white', bordercolor='gray90', borderless=1, focuscolor='', width=150)
    SubModels.grid(row=64, column=1, sticky='W', padx=170, pady=5)
    SubModels.bind('<MouseWheel>', on_vertical)

    lab40 = Label(second_frame,text="  -------------------------------------------------------------------------------------------------------------------------------------------",font="Arial 10", anchor='w', fg='gray70', bg='gray90')
    lab40.grid(row=65, column=1, sticky="W", pady=5)
    lab40.bind('<MouseWheel>', on_vertical)

    subsection3 = Label(second_frame, text="    Graphical settings ", font=("Arial", "16", "italic"), anchor='w', bg='gray90')
    subsection3.grid(row=66, column=1, sticky="W", padx=5, pady=5)
    subsection3.bind('<MouseWheel>', on_vertical)

    lab41 = Label(second_frame, text="* Multiple pages", font=('Arial', 12), anchor='w', bg='gray90')
    lab41.grid(row=67, column=1, sticky="W", padx=5)
    lab41.bind('<MouseWheel>', on_vertical)
    lab42 = Label(second_frame, text="PDF documents with multiple pages. No, Yes ", font="Arial 10", anchor='w', bg='gray90')
    lab42.grid(row=68, column=1, sticky="W", padx=10)
    lab42.bind('<MouseWheel>', on_vertical)
    MultiPages_menu = ttk.Combobox(second_frame, value=MultiPages, state="readonly")
    MultiPages_menu.grid(row=69, column=1, sticky='W', padx=150, pady=5)
    MultiPages_menu.bind('<MouseWheel>', on_vertical)
    lab42_2 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab42_2.grid(row=69, column=1, sticky="W", padx=10, pady=5)
    lab42_2.bind('<MouseWheel>', on_vertical)

    lab43 = Label(second_frame,text="  -------------------------------------------------------------------------------------------------------------------------------------------",font="Arial 10", anchor='w', fg='gray70', bg='gray90')
    lab43.grid(row=70, column=1, sticky="W", pady=5)
    lab43.bind('<MouseWheel>', on_vertical)

    Additional_parameters = ttk.Button(second_frame, text="+ Additional parameters", command=Hide_General_Settings)
    Additional_parameters.grid(row=71, column=1, sticky="W")
    Additional_parameters.bind('<MouseWheel>', on_vertical)
    s.configure('TButton', font=('Arial', 16, 'bold'), background='gray90', relief="flat")
    s.map('TButton', background=[('active', 'gray90')])


    Oplab1 = Label(second_frame, text="* Number of processors", font="Arial 14", anchor='w', bg='gray90')
    Oplab1.bind('<MouseWheel>', on_vertical)
    Oplab2 = Label(second_frame, text="Max number of proccesors selected by default", font="Arial 10", anchor='w', bg='gray90')
    Oplab2.bind('<MouseWheel>', on_vertical)
    Numb_proce = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    Numb_proce.insert(END, '1')
    Numb_proce.bind('<MouseWheel>', on_vertical)
    #balloon = Pmw.Balloon(second_frame)
    #balloon.bind(Oplab1, "Structual simulations and summary stadistics calculations takes a lot of time, \nso is strongly recommended to use the more of processors you can")

    Oplab3 = Label(second_frame, text="Sampling at different times", font=('Arial', 14), anchor='w', bg='gray90')
    Oplab3.bind('<MouseWheel>', on_vertical)
    Oplab4 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    Oplab4.bind('<MouseWheel>', on_vertical)
    sampling = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    sampling.bind('<MouseWheel>', on_vertical)

    Oplab5 = Label(second_frame, text="Generation time", font=('Arial', 14), anchor='w', bg='gray90')
    Oplab5.bind('<MouseWheel>', on_vertical)
    G_time = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    G_time.bind('<MouseWheel>', on_vertical)

    Oplab6 = Label(second_frame, text="Rate of heterogeneity across sites (+G)", font=('Arial', 14), anchor='w',bg='gray90')
    Oplab6.bind('<MouseWheel>', on_vertical)
    Oplab7 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    Oplab7.bind('<MouseWheel>', on_vertical)
    G_menu = ttk.Combobox(second_frame, value=distribution, state="readonly")
    G_menu.config(width=10)
    G_menu.bind("<<ComboboxSelected>>", pick_G)
    G_menu.bind('<MouseWheel>', on_vertical)
    Oplab8 = Label(second_frame, text="      Example:norm                         1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab8.bind('<MouseWheel>', on_vertical)
    G_menu_1 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    G_menu_1.config(width=12)
    G_menu_1.bind('<MouseWheel>', on_vertical)
    Oplab9 = Label(second_frame, text="   (t)", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab9.bind('<MouseWheel>', on_vertical)
    G_menu_2 = ttk.Combobox(second_frame, value=truncated, state="readonly")
    G_menu_2.config(width=3)
    G_menu_2.bind('<MouseWheel>', on_vertical)
    Oplab10 = Label(second_frame, text="1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab10.bind('<MouseWheel>', on_vertical)
    G_menu_3 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    G_menu_3.config(width=12)
    G_menu_3.bind('<MouseWheel>', on_vertical)

    Oplab11 = Label(second_frame, text="Proportion of invariable sites (+I)", font=('Arial', 14), anchor='w',bg='gray90')
    Oplab11.bind('<MouseWheel>', on_vertical)
    Oplab12 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    Oplab12.bind('<MouseWheel>', on_vertical)
    I_menu = ttk.Combobox(second_frame, value=distribution, state="readonly")
    I_menu.config(width=10)
    I_menu.bind("<<ComboboxSelected>>", pick_I)
    I_menu.bind('<MouseWheel>', on_vertical)
    Oplab13 = Label(second_frame, text="      Example:norm                         1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab13.bind('<MouseWheel>', on_vertical)
    I_menu_1 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    I_menu_1.config(width=12)
    I_menu_1.bind('<MouseWheel>', on_vertical)
    Oplab14 = Label(second_frame, text="   (t)", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab14.bind('<MouseWheel>', on_vertical)
    I_menu_2 = ttk.Combobox(second_frame, value=truncated, state="readonly")
    I_menu_2.config(width=3)
    I_menu_2.bind('<MouseWheel>', on_vertical)
    Oplab15 = Label(second_frame, text="1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab15.bind('<MouseWheel>', on_vertical)
    I_menu_3 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    I_menu_3.config(width=12)
    I_menu_3.bind('<MouseWheel>', on_vertical)

    Oplab16 = Label(second_frame, text="Exponential growth rate or Demographic periods", font=('Arial', 14), anchor='w', bg='gray90')
    Oplab16.bind('<MouseWheel>', on_vertical)
    Oplab17 = Label(second_frame, text="Exponential growth rate (0) or demographic periods (1)", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab17.bind('<MouseWheel>', on_vertical)
    GR_menu = ttk.Combobox(second_frame, value=GrowthRate, state="readonly")
    GR_menu.current(0)
    GR_menu.bind("<<ComboboxSelected>>", pick_GR)
    GR_menu.bind('<MouseWheel>', on_vertical)
    Oplab18 = Label(second_frame, text="Exponential growth per individual per generation", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab18.bind('<MouseWheel>', on_vertical)
    Oplab19 = Label(second_frame, text="Example:       1e-5", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab19.bind('<MouseWheel>', on_vertical)
    GR_entry_1 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    GR_entry_1.bind('<MouseWheel>', on_vertical)
    GR_entry_1.configure(state="disabled")
    Oplab20 = Label(second_frame, text="Number periods + N0 + Nend + + duration per period", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab20.bind('<MouseWheel>', on_vertical)
    Oplab21 = Label(second_frame, text="3 1000 1250 1000 1300 1550 2000 1560 1000 3000", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab21.bind('<MouseWheel>', on_vertical)
    GR_entry_2 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    GR_entry_2.bind('<MouseWheel>', on_vertical)
    GR_entry_2.configure(state="disabled")

    Oplab22 = Label(second_frame, text="Migration model and population structure", font=('Arial', 14), anchor='w', bg='gray90')
    Oplab22.bind('<MouseWheel>', on_vertical)
    Oplab23 = Label(second_frame, text="See Manual for more information. Example: 2 2 3 3", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab23.bind('<MouseWheel>', on_vertical)
    MM_entry = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    MM_entry.bind('<MouseWheel>', on_vertical)

    Oplab24 = Label(second_frame, text="Migration rate", font=('Arial', 14), anchor='w', bg='gray90')
    Oplab24.bind('<MouseWheel>', on_vertical)
    Oplab25 = Label(second_frame, text="See Manual for more information. Example: 3 100 800 0.002 0.001 0.003", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab25.bind('<MouseWheel>', on_vertical)
    MR_entry = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    MR_entry.bind('<MouseWheel>', on_vertical)

    Oplab26 = Label(second_frame, text="Events of convergence of demes", font=('Arial', 14), anchor='w', bg='gray90')
    Oplab26.bind('<MouseWheel>', on_vertical)
    Oplab27 = Label(second_frame, text="See Manual for more information. Example: 3 1 2 400 3 4 1900 5 6 2000", font=('Arial', 10), anchor='w', bg='gray90')
    Oplab27.bind('<MouseWheel>', on_vertical)
    CD_entry = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    CD_entry.bind('<MouseWheel>', on_vertical)

    #ConvergenceDemes = ' '.join(ConvergenceDemes)

    lab44 = Label(second_frame,text="  -------------------------------------------------------------------------------------------------------------------------------------------",font="Arial 10", anchor='w', fg='gray70', bg='gray90')
    lab44.grid(row=600, column=1, sticky="W", pady=5)
    lab44.bind('<MouseWheel>', on_vertical)

    save_button = Button(second_frame, text="Save and Continue ➔", command=lambda:save_selected_values(),  bg='white', borderless=1)
    save_button.config(font=('Arial', 14, 'bold'), width=210, height=30)
    save_button.grid(row=700, column=1, pady=5, sticky='W', padx=150)
    save_button.bind('<MouseWheel>', on_vertical)



#Splash Screen timer
splash_root.after(2000, main_window)





mainloop()
