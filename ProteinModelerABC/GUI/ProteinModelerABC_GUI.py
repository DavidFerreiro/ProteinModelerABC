import threading
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
#from tkmacosx import Button
import subprocess as sub
import time
import os
import sys


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)


splash_root = Tk()
file1 = resource_path("Mini.png")
file2 = resource_path("Splash.png")
#img = Image("photo", file="Mini.png")
img = Image("photo", file=file1)
splash_root.call('wm','iconphoto', splash_root._w, img)
#img2 = PhotoImage(file="Splash.png")
img2 = PhotoImage(file=file2)

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
#canvas.master.overrideredirect(True)
canvas.master.title('ProteinModelerABC')
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
Tree_directory = ''
haploid_value = ''
LB = ''

def main_window():
    splash_root.destroy()

    root = Tk()
    root.title('ProteinModelerABC')
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
    global is_oon
    is_oon = True

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
        "Haploid",
        "Diploid"
    ]

    distribution = [
        "fix",
        "uniform",
        "norm",
        "exp",
        "gamma",
        "beta",
        "dirichlet"
    ]

    amino_distribution = [
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
        "rejection",
        "mnlogistic",
        "neuralnet"
    ]

    MultiPages = [
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
        Boton_pop = Button(pop, text="OK", command=pop.destroy, bg="white")
        Boton_pop.grid(row=2, column=2, pady=10)
        
    def popup2(tit, messag):
        pop = Toplevel()
        pop.title(tit)
        #pop.geometry("350x125")
        pop.config(bg="white")
        pop.geometry(f"450x125+{root.winfo_x()+70}+{root.winfo_y()+200}")

        pop_ico = Label(pop, text="🥳", font=('Arial', 40), fg='red', bg='white')
        #pop_ic2 = Label(pop, text="⚠", font=('Arial', 40), fg='red', bg='gray90')
        #pop_ic3 = Label(pop, text="🚨", font=('Arial', 40), fg='red', bg='gray90')
        pop_ico.grid(row=1, column=1, pady=10, padx=10)
        pop_label1 = Label(pop, text = messag, font=('Arial', 14), bg='white')
        pop_label1.grid(row=1, column=2, padx=10)
        Boton_pop = Button(pop, text="OK", command=pop.destroy, bg="white")
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
        newWindows_valores.set("DGREM_Mean DGREM_sd SegSites Grantham_mean_Position Grantham_sd_Position Grantham_sk_Position Grantham_ku_Position")
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

            newWindows.destroy()

        Select = Button(newWindows_frame, text="Select", command=select, bg='gray90', highlightbackground='gray90')
        Select.grid(column=0, row=1, padx=85, pady=5)

    def openEM():
        newWindowss = Toplevel(second_frame)
        newWindowss.title("Empirical Substitution Models")
        #newWindows.geometry("285x315")
        newWindowss.geometry(f"285x160+{root.winfo_x()+110}+{root.winfo_y()}")
        newWindowss.config(background="gray90")
        newWindowss.resizable(False, False)
        newWindowss_frame = ttk.Frame(newWindowss, padding=(6, 3, 12, 12))
        newWindowss_frame.grid(column=0, row=0, sticky=(N, S, E, W))
        #s.configure('TFrame', background='gray90')
        # valores
        newWindowss_valores = StringVar()
        newWindowss_valores.set("Blosum62 CpRev Dayhoff DayhoffDCMUT HIVb HIVw JTT JonesDCMUT LG Mtart Mtrev24 RtRev VT WAG")
        # crear la lista
        newWindowss_lstbox = Listbox(newWindowss_frame, listvariable=newWindowss_valores, selectmode=MULTIPLE, width=30, height=7)
        #newWindowss_lstbox.select_set(first=0, last=7)
        newWindowss_lstbox.grid(column=0, row=0, columnspan=2)
        # frame.bind('<MouseWheel>', no_op)

        # selecionar
        def select():
            global S_model1
            global LB
            S_models = []
            seleccion = newWindowss_lstbox.curselection()
            for i in seleccion:
                entrada = newWindowss_lstbox.get(i)
                S_models.append(entrada)
            S_model1 = str(S_models)[1:-1]
            S_model1 = S_model1.replace("'", "")
            S_model1 = S_model1.replace(",", "")
            LB = S_model1

            if LB != '':
                A_frequencies_menu.configure(state="normal")
                A_frequencies.configure(state="normal")
            else:
                A_frequencies.configure(state="disabled")
                A_frequencies_menu.configure(state="disabled")


            newWindowss.destroy()

        Select = Button(newWindowss_frame, text="Select", command=select, bg='gray90', highlightbackground='gray90')
        Select.grid(column=0, row=1, padx=85, pady=5)

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

    def Name_Template():
        global Tempname
        global Template_directory
        Tempname = Label(second_frame, text=" ", bg="gray90")
        Template_directory = filedialog.askopenfilename()
        if Template_directory == '':
            Tempname.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            Tempname.grid(row=10, column=1, sticky="W", padx=5,)
            lab6_3.config(text='')
            popup("ERROR","You must upload an alignment file")
        else:
            Template_directory = Template_directory.split('/')[-1]
            lab6_3.config(text=Template_directory)
            Tempname.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            Tempname.grid(row=10, column=1, sticky="W", padx=5)

    def Name_Tree():
        global Treename
        global Tree_directory
        Treename = Label(second_frame, text=" ", bg="gray90")
        Tree_directory = filedialog.askopenfilename()
        if Tree_directory == '':
            Treename.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            Treename.grid(row=35, column=1, sticky="W", padx=5)
            lab24_2.config(text='')
            popup("ERROR","You did not upload any file")
        else:
            Tree_directory = Tree_directory.split('/')[-1]
            lab24_2.config(text=Tree_directory)
            Treename.grid(row=35, column=1, sticky="W", padx=5)

    def checkF():
        if vF.get() == 0:
            NPOP.configure(state="disabled")
        else:
            NPOP.configure(state="normal")

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
            SR_menu_2.configure(state="normal")
            SR_menu_3.delete(first=0,last=1000)
            SR_menu_3.configure(state="normal")
        if SR_menu.get() == "exp":
            SR_menu_1.delete(first=0,last=1000)
            SR_menu_2.current(0)
            SR_menu_2.configure(state="normal")
            SR_menu_3.delete(first=0,last=1000)
            SR_menu_3.configure(state="normal")
        if SR_menu.get() == "gamma":
            SR_menu_1.delete(first=0,last=1000)
            SR_menu_2.current(0)
            SR_menu_2.configure(state="normal")
            SR_menu_3.delete(first=0,last=1000)
            SR_menu_3.configure(state="normal")
        if SR_menu.get() == "beta":
            SR_menu_1.delete(first=0,last=1000)
            SR_menu_2.current(0)
            SR_menu_2.configure(state="normal")
            SR_menu_3.delete(first=0,last=1000)
            SR_menu_3.configure(state="normal")
        if SR_menu.get() == "dirichlet":
            SR_menu_1.delete(first=0,last=1000)
            SR_menu_2.current(0)
            SR_menu_2.configure(state="disabled")
            SR_menu_3.delete(first=0,last=1000)
            SR_menu_3.configure(state="disabled")

    A_frequenciesModel = 'fix '
    def pick_amino(e):
        global A_frequenciesModel
        #print('Hola')
        #print(LB)
        #if LB != '':
            #print('adios')
            #A_frequencies_menu.configure(state="normal")
            #A_frequencies.configure(state="normal")
        if A_frequencies_menu.get() == "fix":
            A_frequencies.delete(first=0, last=1000)
            A_frequencies.insert(0, "0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05")
            A_frequenciesModel = 'fix '
        elif A_frequencies_menu.get() == "dirichlet":
            A_frequencies.delete(first=0, last=1000)
            A_frequencies.insert(0, "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1")
            A_frequenciesModel = 'dirichlet '
        #else:
            #A_frequencies.configure(state="disabled")
            #A_frequencies_menu.configure(state="disabled")

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
        if v2.get() == 0:
            lab13.grid(row=35, column=1, sticky="W", padx=5)
            lab13_2.grid(row=36, column=1, sticky="W", padx=5)
            haploid_menu.grid(row=36, column=1, sticky='W', padx=150, pady=5)
            lab14.grid(row=37, column=1, sticky="W", padx=5)
            lab14_2.grid(row=38, column=1, sticky="W", padx=10)
            SR_menu.grid(row=39, column=1, sticky="W", padx=10, pady=5)
            lab14_3.grid(row=38, column=1, sticky="W", padx=10)
            lab14_4.grid(row=38, column=1, sticky="W", padx=10)
            SR_menu_1.grid(row=39, column=1, sticky="W", padx=140, pady=5)
            lab14_5.grid(row=38, column=1, sticky="W", padx=285)
            SR_menu_2.grid(row=39, column=1, sticky="W", padx=285, pady=5)
            lab14_6.grid(row=38, column=1, sticky="W", padx=360)
            SR_menu_3.grid(row=39, column=1, sticky="W", padx=360, pady=5)
            lab15.grid(row=40, column=1, sticky="W", padx=5)
            lab15_2.grid(row=41, column=1, sticky="W", padx=5)
            Ne.grid(row=41, column=1, sticky='W', padx=150, pady=5)
            Additional_parameter.grid(row=42, column=1, sticky="W", padx=5)
            Treename.grid_forget()
            lab24.grid_forget()
            lab24_2.grid_forget()
            upload_PT.grid_forget()
            updateScrollRegion()

        else:
            lab13.grid_forget()
            lab13_2.grid_forget()
            haploid_menu.grid_forget()
            lab14.grid_forget()
            lab14_2.grid_forget()
            lab14_7.grid_forget()
            SR_menu.grid_forget()
            lab14_3.grid_forget()
            lab14_4.grid_forget()
            SR_menu_1.grid_forget()
            lab14_5.grid_forget()
            SR_menu_2.grid_forget()
            lab14_6.grid_forget()
            SR_menu_3.grid_forget()
            lab15.grid_forget()
            lab15_2.grid_forget()
            Ne.grid_forget()
            Additional_parameter.grid_forget()
            Treename.grid(row=35, column=1, sticky="W", padx=5)
            lab24.grid(row=34, column=1, sticky="W", padx=150)
            lab24_2.grid(row=35, column=1, sticky="W", padx=150)
            upload_PT.grid(row=35, column=1, sticky="W", padx=380)
            updateScrollRegion()

    def Hide_SM_Settings():
        global is_oon
        # determin is on ot off
        if is_oon:
            Additional_parameters.config(text="- Optional parameters")
            is_oon = False
            lab28.grid(row=82, column=1, sticky='W', padx=5)
            lab28_2.grid(row=83, column=1, sticky='W', padx=10)
            lab28_3.grid(row=84, column=1, sticky='W', padx=10)
            G_menu.grid(row=85, column=1, sticky='W', padx=10, pady=5)
            G_menu_1.grid(row=85, column=1, sticky='W', padx=140, pady=5)
            lab28_4.grid(row=84, column=1, sticky='W', padx=285)
            G_menu_2.grid(row=85, column=1, sticky='W', padx=285, pady=5)
            lab28_5.grid(row=84, column=1, sticky='W', padx=360)
            G_menu_3.grid(row=85, column=1, sticky='W', padx=360, pady=5)
            lab29.grid(row=86, column=1, sticky='W', padx=5)
            lab29_2.grid(row=87, column=1, sticky='W', padx=10)
            I_menu.grid(row=89, column=1, sticky='W', padx=10, pady=5)
            lab29_3.grid(row=88, column=1, sticky='W', padx=10)
            I_menu_1.grid(row=89, column=1, sticky='W', padx=140, pady=5)
            lab29_4.grid(row=88, column=1, sticky='W', padx=285)
            I_menu_2.grid(row=89, column=1, sticky='W', padx=285, pady=5)
            lab29_5.grid(row=88, column=1, sticky='W', padx=360)
            I_menu_3.grid(row=89, column=1, sticky='W', padx=360, pady=5)
            updateScrollRegion()
        else:
            Additional_parameters.config(text="+ Optional parameters")
            is_oon = True
            lab28.grid_forget()
            lab28_2.grid_forget()
            G_menu.grid_forget()
            lab28_3.grid_forget()
            G_menu_1.grid_forget()
            lab28_4.grid_forget()
            G_menu_2.grid_forget()
            lab28_5.grid_forget()
            G_menu_3.grid_forget()
            lab29.grid_forget()
            lab29_2.grid_forget()
            I_menu.grid_forget()
            lab29_3.grid_forget()
            I_menu_1.grid_forget()
            lab29_4.grid_forget()
            I_menu_2.grid_forget()
            lab29_5.grid_forget()
            I_menu_3.grid_forget()
            updateScrollRegion()

    def Hide_col_Settings():
        global is_on
        # determin is on ot off
        if is_on:
            Additional_parameter.config(text="- Optional parameters")
            is_on = False
            lab16.grid(row=43, column=1, sticky="W", padx=5)
            lab16_2.grid(row=44, column=1, sticky="W", padx=10)
            sampling.grid(row=45, column=1, sticky="W", padx=150)
            lab17.grid(row=46, column=1, sticky="W", padx=5)
            lab17_2.grid(row=47, column=1, sticky="W", padx=10)
            G_time.grid(row=48, column=1, sticky="W", padx=150)
            lab18.grid(row=49, column=1, sticky="W", padx=5)
            lab18_2.grid(row=50, column=1, sticky="W", padx=10)
            GR_menu.grid(row=51, column=1, sticky="W", padx=150, pady=5)
            lab19.grid(row=52, column=1, sticky="W", padx=5)
            lab19_2.grid(row=53, column=1, sticky="W", padx=10)
            GR_entry_1.grid(row=54, column=1, sticky="W", padx=30, pady=5)
            lab20.grid(row=52, column=1, sticky="W", padx=240)
            lab20_2.grid(row=53, column=1, sticky="W", padx=240)
            GR_entry_2.grid(row=54, column=1, sticky="W", padx=250, pady=5)
            lab21.grid(row=55, column=1, sticky="W", padx=5)
            lab21_2.grid(row=56, column=1, sticky="W", padx=10)
            MM_entry.grid(row=57, column=1, sticky="W", padx=150, pady=5)
            lab22.grid(row=58, column=1, sticky="W", padx=5)
            lab22_2.grid(row=59, column=1, sticky="W", padx=10)
            MR_entry.grid(row=60, column=1, sticky="W", padx=150, pady=5)
            lab23.grid(row=61, column=1, sticky="W", padx=5)
            lab23_2.grid(row=62, column=1, sticky="W", padx=10)
            CD_entry.grid(row=63, column=1, sticky="W", padx=150, pady=5)
            updateScrollRegion()
        else:
            Additional_parameter.config(text="+ Optional parameters")
            is_on = True
            lab16.grid_forget()
            lab16_2.grid_forget()
            sampling.grid_forget()
            lab17.grid_forget()
            lab17_2.grid_forget()
            G_time.grid_forget()
            lab18.grid_forget()
            lab18_2.grid_forget()
            GR_menu.grid_forget()
            lab19.grid_forget()
            lab19_2.grid_forget()
            GR_entry_1.grid_forget()
            lab20.grid_forget()
            lab20_2.grid_forget()
            GR_entry_2.grid_forget()
            lab21.grid_forget()
            lab21_2.grid_forget()
            MM_entry.grid_forget()
            lab22.grid_forget()
            lab22_2.grid_forget()
            MR_entry.grid_forget()
            lab23.grid_forget()
            lab23_2.grid_forget()
            CD_entry.grid_forget()
            updateScrollRegion()

    def save_selected_values():
        ## Chain
        if ChainEntry.get() != '':
            Chainlab_3.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            Chainlab_3.grid(row=13, column=1, sticky="W", padx=5)
        else:
            Chainlab_3.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            Chainlab_3.grid(row=13, column=1, sticky="W", padx=5)
            popup("ERROR", "Error in number of simulations\nYou must enter the number of simulations")
            sys.tracebacklimit = 0
            raise ValueError()


        ## Number of simulations check
        Numb_simu2 = Numb_simu.get()
        if Numb_simu.get() != '':
            if Numb_simu2.isnumeric() == True:
                if int(Numb_simu2) < 100000000:
                    if int(Numb_simu2) > 0:
                        lab8_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                        lab8_2.grid(row=19, column=1, sticky="W", padx=5)
                    else:
                        lab8_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                        lab8_2.grid(row=19, column=1, sticky="W", padx=5)
                        popup("ERROR", "Error in number of simulations\nYou must enter a value between 0-100000")
                        sys.tracebacklimit = 0
                        raise ValueError()
                else:
                    lab8_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab8_2.grid(row=19, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error in number of simulations\nYou must enter a value between 0-100000")
                    sys.tracebacklimit = 0
                    raise ValueError()
            else:
                lab8_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab8_2.grid(row=19, column=1, sticky="W", padx=5)
                popup("ERROR", "Error in number of simulations\nYou must enter a number")
                sys.tracebacklimit = 0
                raise ValueError()
        else:
            lab8_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab8_2.grid(row=19, column=1, sticky="W", padx=5)
            popup("ERROR", "Error in number of simulations\nYou must enter the number of simulations")
            sys.tracebacklimit = 0
            raise ValueError()

        #Number of processors   Oplab3.grid(row=22, column=1, sticky="W", padx=10)
        Numb_proce2 = Numb_proce.get()
        if Numb_proce.get() != '':
            if Numb_proce2.isnumeric() == True:
                if int(Numb_proce.get()) >= int(os.cpu_count()):
                    lab13_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                    lab13_2.grid(row=22, column=1, sticky="W", padx=5)
                else:
                    lab13_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab13_2.grid(row=22, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error in number of processors\nYou selected more processors than your machine has")
                    sys.tracebacklimit = 0
                    raise ValueError()
            else:
                lab13_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab13_2.grid(row=22, column=1, sticky="W", padx=5)
                popup("ERROR", "Error in number of processors\nProcessors value has to be a integer")
                sys.tracebacklimit = 0
                raise ValueError()
        else:
            lab13_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab13_2.grid(row=22, column=1, sticky="W", padx=5)
            popup("ERROR", "Error in number of processors\nYou must select a number of processors")
            sys.tracebacklimit = 0
            raise ValueError()


        ## COALESCENT
        if v2.get() == 0:
            ## Haploid/Diploid check
            if str(haploid_menu.get()) != '':
                lab13_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                lab13_2.grid(row=36, column=1, sticky="W", padx=9)
            else:
                lab13_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab13_2.grid(row=36, column=1, sticky="W", padx=9)
                popup("ERROR", "Error in Haploid/Diploid simulated data\nYou must select the haploid or diploid option")
                sys.tracebacklimit = 0
                raise ValueError()


            ## Aminoacid substitution rate per site
            if str(SR_menu.get()) == '':
                lab14_7.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab14_7.grid(row=39, column=1, sticky="W", padx=5)
                popup("ERROR", "Error in Aminoacid Substitution Rate\nYou must select at least one model")
                sys.tracebacklimit = 0
                raise ValueError()
            else:
                if str(SR_menu.get()) == 'fix' and str(SR_menu_1.get()) == '':
                    lab14_7.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab14_7.grid(row=39, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error in Aminoacid Substitution Rate\nYou must enter a value")
                    sys.tracebacklimit = 0
                    raise ValueError()
                elif str(SR_menu.get()) == 'norm' and str(SR_menu_2.get()) == 't' and str(SR_menu_3.get()) =='':
                    lab14_7.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab14_7.grid(row=39, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error in Aminoacid Substitution Rate\nIf you select truncated option you must enter a value")
                    sys.tracebacklimit = 0
                    raise ValueError()
                elif str(SR_menu.get()) == 'exp' and str(SR_menu_2.get()) == 't' and str(SR_menu_3.get()) =='':
                    lab14_7.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab14_7.grid(row=39, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error in Aminoacid Substitution Rate\nIf you select truncated option you must enter a value")
                    sys.tracebacklimit = 0
                    raise ValueError()
                elif str(SR_menu.get()) == 'gamma' and str(SR_menu_2.get()) == 't' and str(SR_menu_3.get()) =='':
                    lab14_7.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab14_7.grid(row=39, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error in Aminoacid Substitution Rate\nIf you select truncated option you must enter a value")
                    sys.tracebacklimit = 0
                    raise ValueError()
                elif str(SR_menu.get()) == 'beta' and str(SR_menu_2.get()) == 't' and str(SR_menu_3.get()) =='':
                    lab14_7.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab14_7.grid(row=39, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error in Aminoacid Substitution Rate\nIf you select truncated option you must enter a value")
                    sys.tracebacklimit = 0
                    raise ValueError()
                else:
                    lab14_7.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                    lab14_7.grid(row=39, column=1, sticky="W", padx=5)

            ## Efective population size check
            Ne2 = Ne.get()
            if Ne.get() != '':
                if Ne2.isnumeric() == True:
                    lab15_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                    lab15_2.grid(row=41, column=1, sticky="W", padx=5)
                else:
                    lab15_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab15_2.grid(row=41, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error in effective population size\nYou must enter a integer")
                    sys.tracebacklimit = 0
                    raise ValueError()
            else:
                lab15_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                lab15_2.grid(row=41, column=1, sticky="W", padx=5)
                popup("ERROR", "Error in effective population size\nYou must enter the number of simulations")
                sys.tracebacklimit = 0
                raise ValueError()


        # Substitution Model
        global SS_model
        SS_model = []
        if vF.get() == 1:
            SS_model.append('Fitness')
        if vN.get() == 1:
            SS_model.append('Neutral')
        if vF.get() == 0 and vN.get() == 0:
            lab26_3.config(text="SCS models ❌", bg="gray90")
            popup("ERROR", "Error in effective population size\nYou must enter the number of simulations")
            sys.tracebacklimit = 0
            raise ValueError()
        else:
            lab26_3.config( text="SCS models ✅", bg="gray90")

        SS_model1 = str(SS_model)[1:-1]
        SS_model1 = SS_model1.replace("'", "")
        SS_model1 = SS_model1.replace(",", "")

        #Temp
        TEMP = Temp.get()
        if Temp.get() != '':
            if TEMP.isnumeric() == True:
                lab127_1.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            elif TEMP.replace('.', '', 1).isdigit() == True:
                lab127_1.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            else:
                lab127_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                popup("ERROR", "Error in termodynamic temperature\nYou must enter a float number")
                sys.tracebacklimit = 0
                raise ValueError()
        else:
            lab127_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            popup("ERROR", "Error in termodynamic temperature\nYou must enter a float number")
            sys.tracebacklimit = 0
            raise ValueError()

        #S0
        s0 = S0.get()
        if S0.get() != '':
            if s0.isnumeric() == True:
                lab128_1.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            elif s0.replace('.', '', 1).isdigit() == True:
                lab128_1.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            else:
                lab128_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                popup("ERROR", "Error in unfolded protein \nconfigurational entropy (S0)\nYou must enter a float number")
                sys.tracebacklimit = 0
                raise ValueError()
        else:
            lab128_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            popup("ERROR", "Error in unfolded protein \nconfigurational entropy (S0)\nYou must enter a float number")
            sys.tracebacklimit = 0
            raise ValueError()

        #SC1
        sC1 = SC1.get()
        if SC1.get() != '':
            if sC1.isnumeric() == True:
                lab129_1.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            elif sC1.replace('.', '', 1).isdigit() == True:
                lab129_1.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            else:
                lab129_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                popup("ERROR", "Error in misfolded protein \nconfigurational entropy (SC1)\nYou must enter a float number")
                sys.tracebacklimit = 0
                raise ValueError()
        else:
            lab129_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            popup("ERROR", "Error in misfolded protein \nconfigurational entropy (SC1)\nYou must enter a float number")
            sys.tracebacklimit = 0
            raise ValueError()


        #SC0
        sC0 = SC0.get()
        if SC0.get() != '':
            if sC0.isnumeric() == True:
                lab130_1.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            elif s0.replace('.', '', 1).isdigit() == True:
                lab130_1.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            else:
                lab130_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                popup("ERROR", "Error in misfolded protein \nconfigurational entropy offset (SC0)\nYou must enter a float number")
                sys.tracebacklimit = 0
                raise ValueError()
        else:
            lab130_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            popup("ERROR", "Error in misfolded protein \nconfigurational entropy offset (SC0)\nYou must enter a float number")
            sys.tracebacklimit = 0
            raise ValueError()

        #REM3
        rem3 = REM3.get()
        if REM3.get() != '':
            if rem3.isnumeric() == True:
                lab131_1.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            else:
                lab131_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                popup("ERROR", "Error in the third cumulant in REM \ncalculation to compute the structural \nprotein energy (REM3)\nYou must enter a integer")
                sys.tracebacklimit = 0
                raise ValueError()
        else:
            lab131_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            popup("ERROR", "Error in the third cumulant in REM \ncalculation to compute the structural \nprotein energy (REM3)\nYou must enter a integer")
            sys.tracebacklimit = 0
            raise ValueError()

        #NPOP
        npop = NPOP.get()
        if NPOP.get() != '':
            if npop.isnumeric() == True:
                lab132_1.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            else:
                lab132_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                popup("ERROR", "Error in the population size to simulate \nmolecular evolution under Fitness \nsite-dependent SCS model\nYou must enter a integer")
                sys.tracebacklimit = 0
                raise ValueError()
        else:
            lab132_1.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            popup("ERROR", "Error in the population size to simulate \nmolecular evolution under Fitness \nsite-dependent SCS model\nYou must enter a integer")
            sys.tracebacklimit = 0
            raise ValueError()

        ## ABC iteractions
        ABC_itera2 = ABC_itera.get()
        if ABC_itera2 != "":
            if ABC_itera2.isnumeric() == True:
                if int(ABC_itera2) <= int(Numb_simu2):
                    if int(ABC_itera2) > 0:
                        lab30_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                        lab30_2.grid(row=94, column=1, sticky="W", padx=5)
                    else:
                        lab30_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                        lab30_2.grid(row=94, column=1, sticky="W", padx=5)
                        popup("ERROR","Error ABC iterations\nYou must enter a value greater than 0")
                        sys.tracebacklimit = 0
                        raise ValueError()
                else:
                    lab30_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab30_2.grid(row=94, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error ABC iterations\nYou must enter a value less or equal\n to the number of simulations")
                    sys.tracebacklimit = 0
                    raise ValueError()
            else:
                if ABC_itera2.lstrip("-").isnumeric() == True:
                    lab30_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab30_2.grid(row=94, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error ABC iterations\nYou must enter a positive value")
                    sys.tracebacklimit = 0
                    raise ValueError()
                elif ABC_itera2.replace('.', '', 1).isdigit() == True:
                    lab30_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab30_2.grid(row=94, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error ABC iterations\nYou must enter an integer")
                    sys.tracebacklimit = 0
                    raise ValueError()
                else:
                    lab30_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab30_2.grid(row=94, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error ABC iterations\nYou must enter a number")
                    sys.tracebacklimit = 0
                    raise ValueError()
        else:
            lab27_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab27_2.grid(row=94, column=1, sticky="W", padx=5)
            popup("ERROR", "Error ABC iterations\nYou must enter a value")
            sys.tracebacklimit = 0
            raise ValueError()

        ## ABC tolerance
        ABC_tolerance2 = ABC_tolerance.get()
        if ABC_tolerance2 != "":
            if ABC_tolerance2.isnumeric() == True:
                if float(ABC_tolerance2) <= 1:
                    if float(ABC_tolerance2) > 0:
                        lab31_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                        lab31_2.grid(row=97, column=1, sticky="W", padx=5)
                    else:
                        lab31_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                        lab31_2.grid(row=97, column=1, sticky="W", padx=5)
                        popup("ERROR","Error ABC tolerance\nYou must enter a value greater than 0")
                        sys.tracebacklimit = 0
                        raise ValueError()
                else:
                    lab31_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab31_2.grid(row=97, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error ABC tolerance\nYou must enter a value less than 1")
                    sys.tracebacklimit = 0
                    raise ValueError()
            else:
                if ABC_tolerance2.lstrip("-").isnumeric() == True:
                    lab31_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab31_2.grid(row=97, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error ABC tolerance\nYou must enter a positive value")
                    sys.tracebacklimit = 0
                    raise ValueError()
                elif ABC_tolerance2.replace('.', '', 1).isdigit() == True:
                    lab31_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
                    lab31_2.grid(row=97, column=1, sticky="W", padx=5)
                else:
                    lab31_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
                    lab31_2.grid(row=97, column=1, sticky="W", padx=5)
                    popup("ERROR", "Error ABC tolerance\nYou must enter a number")
                    sys.tracebacklimit = 0
                    raise ValueError()

        else:
            lab29_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab29_2.grid(row=97, column=1, sticky="W", padx=5)
            popup("ERROR", "Error ABC tolerance\nYou must enter a value")
            sys.tracebacklimit = 0
            raise ValueError()

        #ABC Method
        if str(ABCMethod_menu.get()) != '':
            lab32_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            lab32_2.grid(row=100, column=1, sticky="W", padx=5)
        else:
            lab32_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab32_2.grid(row=100, column=1, sticky="W", padx=5)
            popup("ERROR", "Error in ABC method\nYou must select ABC method option")
            sys.tracebacklimit = 0
            raise ValueError()

        ## Summary stadistics
        if SStadistics != "":
            lab33_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            lab33_2.grid(row=103, column=1, sticky="W", padx=9)
        else:
            lab33_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab33_2.grid(row=103, column=1, sticky="W", padx=9)
            popup("ERROR", "Error in Summary Stadistics\nYou must select at least 1 summary stadistics")
            sys.tracebacklimit = 0
            raise ValueError()

        ## Multipages
        if str(MultiPages_menu.get()) != "":
            lab34_2.config(text="✅", font=('Arial', 12), fg='green', bg='gray90')
            lab34_2.grid(row=106, column=1, sticky="W", padx=9)
        else:
            lab34_2.config(text="❌", font=('Arial', 12), fg='red', bg='gray90')
            lab34_2.grid(row=106, column=1, sticky="W", padx=9)
            popup("ERROR", "Error in multiple pages\nYou must select a multiple pages value")
            sys.tracebacklimit = 0
            raise ValueError()


        file = open('Settings.txt', 'w')
        file.write('######################################################################################\n')
        file.write('#####                                                                             ####\n')
        file.write('#####   Settings file for ProteinModelerABC                                       ####\n')
        file.write('#####   Selection of the best-fitting empirical or structural substitution model  ####\n')
        file.write('#####   of protein evolution by approximate Bayesian computation                  ####\n')
        file.write('#####   David Ferreiro, Catarina Branco, Ugo Bastolla and Miguel Arenas           ####\n')
        file.write('#####   (c) 2023                                                                  ####\n')
        file.write('#####   Contact: david.ferreiro.garcia@uvigo.es / ferreirogarciadavid@gmail.com	  ####\n')
        file.write('#####	                                                                          ####\n')
        file.write('##### 	x-------------x-------------x-------------x-------------x-------------x   ####\n')
        file.write('#####   Parameters with an "*" are mandatory (need to be specified)               ####\n')
        file.write('#####   Parameter values must be introduced immediately after the "="             ####\n')
        file.write('##### 	x-------------x-------------x-------------x-------------x-------------x   ####\n')
        file.write('#####                                                                             ####\n')
        file.write('######################################################################################\n\n\n')
        file.write('########################################\n')
        file.write('## GENERAL INPUT DATA AND INFORMATION ##\n')
        file.write('########################################\n\n')
        file.write('# File with the target alignment of protein sequences. Phylip format, see documentation for details -- MANDATORY --\n')
        file.write('*NameOfPhylipFile=' + str(alignment_directory2) + '\n\n')
        file.write('# Consideration of indels. "Ignored" (indels are ignored), "NewState" (indels are considered as a new state) -- MANDATORY --\n')
        file.write('*Indels=' + indells_menu.get() + '\n\n')
        file.write('# File with the protein structure (PDB file). Protein structure used for structurally constrained substitution models and certain summary statistics -- MANDATORY --\n')
        file.write('*Template=' + str(Template_directory) + '\n\n')
        file.write('# Chain of the protein structure. See documentation for details -- MANDATORY --\n')
        file.write('*Chain=' + ChainEntry.get() + '\n\n')
        file.write('# Show running information (simulations and summary statistics) on the screen. It increases the computer time. See documentation for details -- MANDATORY --\n')
        file.write('*ShowInformationScreen=' + Show_RI_menu.get() + '\n\n\n')
        file.write('#######################################\n')
        file.write('## SETTINGS FOR THE SIMULATION PHASE ##\n')
        file.write('#######################################\n\n')
        file.write('                           GENERAL SIMULATION SETTINGS\n')
        file.write('------------------------------------------------------------------------------------\n\n')
        file.write('# Total number of simulations -- MANDATORY --\n')
        file.write('*NumberOfSimulations=' + Numb_simu.get() + '\n\n')
        file.write('# Number of available processors to run the simulations in parallel (1 in case of running without parallelization). See documentation for details -- MANDATORY --\n')
        file.write('*NumberOfProcessors=' + Numb_proce.get() + '\n\n')
        file.write('# Save simulated data. We recommended do not save the simulated data because it requires a lot of space -- MANDATORY --\n')
        file.write('*SaveSimulations=' + Save_S_menu.get() + '\n' + '\n')
        file.write('------------------------------------------------------------------------------------\n')
        file.write('                                EVOLUTIONARY HISTORY\n')
        file.write('------------------------------------------------------------------------------------\n')
        file.write('### The user should select a coalescent simulation (with user-specified parameters) or a rooted phylogenetic tree (provided by the user) upon which protein evolution is simulated ###\n\n')
        file.write('# Simulate the evolutionary histories with the coalescent “Coal” or a phylogenetic tree is specified “Phylo”. See documentation for details -- MANDATORY --\n')
        if v2.get() == 0:
            coal = 'Coal'
        elif v2.get() == 1:
            coal = 'Phylo'
        file.write('*CoalescentOrPhylogeny=' + str(coal) + '\n\n')
        file.write(' 				Coalescent parameters\n')
        file.write('|"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""\n')
        file.write('|# Haploid or Diploid data. Haploid=1, Diploid=2 -- MANDATORY IF COALESCENT --\n')
        if haploid_menu.get() == 'Haploid':
            haploid_value = 1
        elif haploid_menu.get() == 'Diploid':
            haploid_value = 2
        else:
            haploid_value=''
        file.write('|*Haploid/Diploid=' + str(haploid_value) + '\n')
        file.write('|\n')
        file.write('|# Amino acid substitution rate. It can be fixed “fix” or include a prior distribution: uniform, gamma, beta, normal, exponential. See documentation for details -- MANDATORY IF COALESCENT --\n')
        file.write('|*SubstitutionRate=' + SR_menu.get() + ' ' + SR_menu_1.get() + ' ' + SR_menu_2.get() + SR_menu_3.get() + '\n')
        file.write('|\n')
        file.write('|# Population size -- MANDATORY IF COALESCENT --\n')
        file.write('|*PopulationSize=' + Ne.get() + '\n')
        file.write('|\n')
        file.write('|# Longitudinal sampling. Requires GenerationTime. See documentation for details\n')
        file.write('|DatedTips=' + sampling.get() + '\n')
        file.write('|\n')
        file.write('|# Generation time. See documentation for details\n')
        file.write('|GenerationTime=' + G_time.get() + '\n')
        file.write('|\n')
        file.write('|# Exponential growth rate or demographic periods. See documentation for details\n')
        file.write('|GrowthRate=' + str(GR_menu.get()) + ' ' + str(GR_entry_1.get()) + str(GR_entry_2.get()) + '\n')
        file.write('|\n')
        file.write('|# Migration model and population structure. See documentation for details\n')
        file.write('|MigrationModel=' + MM_entry.get() + '\n')
        file.write('|\n')
        file.write('|# Migration rate. See documentation for details\n')
        file.write('|MigrationRate=' + MR_entry.get() + '\n')
        file.write('|\n')
        file.write('|# Events of convergence of demes. See documentation for details\n')
        file.write('|ConvergenceDemes=' + CD_entry.get() + '\n')
        file.write('|\n')
        file.write('| 			     Rooted phylogenetic tree\n')
        file.write('|"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""\n')
        file.write('|# File with user-specified phylogenetic tree. Newark format. See documentation for details -- MANDATORY IF USER-SPECIFIED PHYLOGENETIC TREE --\n')
        file.write('|*Tree=' + str(Tree_directory) + '\n')
        file.write('|\n')
        file.write('------------------------------------------------------------------------------------\n')
        file.write('				SUBSTITUTION MODEL\n')
        file.write('------------------------------------------------------------------------------------\n')
        file.write('# Substitution models of protein evolution that are evaluated. Select at least one structurally constrained substitution (SCS) model (Fitness or Neutral) and one desired empirical substitution model that ideally should be the best-fitting substitution model selected with ProtTest or other framework (i.e., Blosum62, CpRev, Dayhoff, DayhoffDCMUT, HIVb, HIVw, JTT, JonesDCMUT, LG, Mtart, Mtmam, Mtrev24, RtRev, VT, WAG). The empirical models should be specified before the SCS models if any empirical model is selected. See documentation for details -- MANDATORY --\n')
        file.write('*SubstitutionModel=' + str(S_model1) + ' ' + str(SS_model1) + '\n\n')
        file.write(' 			   Empirical substitution model\n')
        file.write('|"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""\n')
        file.write('|# Amino acid frequencies. See documentation for details -- MANDATORY IF EMPIRICAL MODEL --\n')
        file.write('|*AminoacidFrequencies=' + str(A_frequenciesModel) + A_frequencies.get() + '\n')
        file.write('|\n')
        file.write('|# Rate of heterogeneity across sites, +G. It can be fixed “fix” or include a prior distribution: uniform, gamma, beta, normal, exponential\n')
        file.write('|RateHetSites=' + G_menu.get() + ' ' + G_menu_1.get() + ' ' + G_menu_2.get() + G_menu_3.get() + '\n')
        file.write('|\n')
        file.write('|# Proportion of invariable sites, +I. It can be fixed “fix” or include a prior distribution: uniform, gamma, beta, normal, exponential\n')
        file.write('|PropInvSites=' + I_menu.get() + ' ' + I_menu_1.get() + ' ' + I_menu_2.get() + I_menu_3.get() + '\n')
        file.write('|\n')
        file.write('| 			          SCS model/s\n')
        file.write('|"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""\n')
        file.write('|# Thermodynamic temperature -- MANDATORY IF SCS MODEL --\n')
        file.write('|*TEMP=' + Temp.get() + '\n')
        file.write('|\n')
        file.write('|# Configurational entropies per residue (unfolded) and offset (misfolded). See documentation for details -- MANDATORY IF SCS MODEL --\n')
        file.write('|*S0=' + S0.get()+ '\n')
        file.write('|*SC1=' + SC1.get() + '\n')
        file.write('|*SC0=' + SC0.get() + '\n')
        file.write('|\n')
        file.write('|# Third cumulant in REM calculation -- MANDATORY IF SCS MODEL --\n')
        file.write('|*REM3=' + REM3.get() + '\n')
        file.write('|\n')
        file.write('|# Population size considered for the SCS model -- MANDATORY IF FITNESS SCS MODEL -- \n')
        file.write('|*NPOP=' + NPOP.get() + '\n')
        file.write('|\n')
        file.write('#######################################\n')
        file.write('## SETTINGS FOR THE ESTIMATION PHASE ##\n')
        file.write('#######################################\n\n')
        file.write('# ABC iterations used for the cross validation. See documentation for details -- MANDATORY --\n')
        file.write('*ABCIterations=' + ABC_itera2 + '\n\n')
        file.write('# ABC tolerance. Proportion of simulations closest to real data to retain in the ABC procedure.  See documentation for details -- MANDATORY --\n')
        file.write('*ABCTolerance=' + ABC_tolerance2 + '\n\n')
        file.write('# ABC method (rejection, mnlogistic or neuralnet).  See documentation for details -- MANDATORY -- \n')
        file.write('*ABCMethod=' + str(ABCMethod_menu.get()) + '\n\n')
        file.write('# Summary statistics that are used for the ABC estimation. See documentation for details -- MANDATORY --\n')
        file.write('*SummaryStatistics=' + SStadistics + '\n\n')
        file.write('# PDF documents with multiple output plots per page -- MANDATORY --\n')
        file.write('*MultiPage=' + str(MultiPages_menu.get()) + '\n\n')
        file.close()



        # Toplevel object which will be treated as a new window
        Running_window = Toplevel(second_frame)
        # sets the title of the Toplevel widget
        Running_window.title("Running ProteinModelerABC")
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
        lab_r_3.place(x=24, y=80)

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
                #popup("DONE!!\n, ProteinModelerABC has finished correctly ")
                
        def reanalyze():
            pass

        #os.system('chmod +x ProteinModelerABCGeneral.py')
        text = Text(Running_window, height=25, width=60, highlightthickness=1, highlightbackground='black',borderwidth=10)
        text.grid(row=3, column=1, columnspan=2, padx=15, pady=35)
        
        def actualizartext():
            threading.Timer(1.0, actualizartext).start()
            
        def start_ProteinModelerABC():
            #print(os.getcwd())
            Progress['value'] = 0
            lab_r_2.configure(text="0%")
            call()
            Progress.start(20)
            os.system('chmod +x ProteinModelerABC.py')
            #threading.Timer(1.0, start_ProteinModelerABC).start()
            #global p
            #p = sub.Popen('./ProteinModelerABC-1.py', stdout=sub.PIPE, stderr=sub.PIPE, bufsize=0) #sin stderr
            p = sub.Popen('./ProteinModelerABC.py', stdout=sub.PIPE, bufsize=0, universal_newlines=True) #bufsize=1 sin stderr
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
                popup("STOPED\n", "ProteinModelerABC has been stopped")
                lab_r_2.configure(text="0%")
                Progress.config(mode='determinate')
                Progress['value']=0
                stopbutton.destroy()
                

            stopbutton = Button(Running_window, text="Stop", font=('Avenir Next', 12, 'bold'), command=stopsubprocess, bg='white', fg='black', highlightbackground='gray90')  # command=start, relief='sunken', bg='white', bordercolor='darkblue', borderless=1, focuscolor='', width=50)
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
                        popup2("DONE!!\n", "ProteinModelerABC has finished correctly")
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

        # p = sub.Popen('./ProteinModelerABC-2.py',stdout=sub.PIPE,stderr=sub.PIPE)
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

        # p = sub.Popen('./ProteinModelerABC-3.py',stdout=sub.PIPE,stderr=sub.PIPE)
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

        # p = sub.Popen('./ProteinModelerABC-4.py',stdout=sub.PIPE,stderr=sub.PIPE)
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

        start = Button(Running_window, text="Start", font=('Avenir Next', 12, 'bold'), command=lambda: threading.Thread(target=start_ProteinModelerABC).start(), relief='flat', bg='lightblue', fg='black', highlightbackground='gray90')  # command=start, relief='sunken', bg='white', bordercolor='darkblue', borderless=1, focuscolor='', width=50)
        start.place(x=390, y=80)
        #start.focus()

        Running_window.mainloop()


    #Pmw.initialise(root)
    sep0 = Label(second_frame, text="  -------------------------------------------------------------------------------------------------------------------------------------------",font="Arial 10", anchor='w', fg='gray70', bg='gray90')
    sep0.grid(row=0, column=1, sticky="W")
    sep0.bind('<MouseWheel>', on_vertical)

    lab2 = Label(second_frame, text=" ## GENERAL INPUT DATA AND INFORMATION ## ", font=('Arial', 16, 'bold'),bg='gray90')
    lab2.grid(row=1, column=1, sticky="W", padx=5, pady=5)
    lab2.bind('<MouseWheel>', on_vertical)

    #subsection1 = Label(second_frame, text="   ### Settings for the simulation phase ###", font=("Arial", "14", "italic"), anchor='w', bg='gray90')
    #subsection1.grid(row=2, column=1, sticky="W", padx=5, pady=5)
    #subsection1.bind('<MouseWheel>', on_vertical)

    lab3 = Label(second_frame, text="* Name of input file with the protein MSA:", font="Arial 14", bg='gray90')
    lab3.grid(row=3, column=1, sticky="W", padx=5)
    lab3.bind('<MouseWheel>', on_vertical)
    lab3_2 = Label(second_frame, text="", font="Arial 14", anchor='w', bg='gray90')
    lab3_2.grid(row=4, column=1, sticky="W", padx=35)
    lab3_2.bind('<MouseWheel>', on_vertical)
    upload = Button(second_frame, text="Upload ↑", command=Name_file, bg='white',highlightbackground='gray90')
    upload.grid(row=4, column=1, sticky="W", padx=380)
    upload.bind('<MouseWheel>', on_vertical)

    lab4 = Label(second_frame, text="* How to consider indels (gaps):", font="Arial 14", anchor='w', bg='gray90')
    lab4.grid(row=5, column=1, sticky="W", padx=5)
    lab4.bind('<MouseWheel>', on_vertical)
    lab5 = Label(second_frame,text="Ignored (by default), New state (indels are considered as a new state):",font="Arial 10", anchor='w', bg='gray90')
    lab5.grid(row=6, column=1, sticky="W", padx=10)
    lab5.bind('<MouseWheel>', on_vertical)
    indells_menu = ttk.Combobox(second_frame, value=gaps, state="readonly")
    indells_menu.current(0)
    indells_menu.grid(row=7, column=1, sticky='W', padx=150,pady=5)
    indells_menu.bind('<MouseWheel>', on_vertical)

    lab6 = Label(second_frame, text="* Representative protein structure (PDB file):", anchor='w', font="Arial 14", bg='gray90')
    lab6.grid(row=8, column=1, sticky="W", padx=5)
    lab6.bind('<MouseWheel>', on_vertical)
    lab6_2 = Label(second_frame,text="Used in the SCS models simulations and to predict proteins free energy",font="Arial 10", anchor='w', bg='gray90')
    lab6_2.grid(row=9, column=1, sticky="W", padx=10)
    lab6_2.bind('<MouseWheel>', on_vertical)
    lab6_3 = Label(second_frame, text="", anchor='w', font="Arial 14", bg='gray90')
    lab6_3.grid(row=10, column=1, sticky="W", padx=35)
    lab6_3.bind('<MouseWheel>', on_vertical)
    upload_T = Button(second_frame, text="Upload ↑", command=Name_Template, bg='white', highlightbackground='gray90')
    upload_T.grid(row=10, column=1, sticky="W", padx=380)
    upload_T.bind('<MouseWheel>', on_vertical)

    Chainlab = Label(second_frame, text="* Chain of the representative PDB:", anchor='w', font="Arial 14", bg="gray90")
    Chainlab.grid(row=11, column=1,  sticky="W", padx=5)
    Chainlab.bind('<MouseWheel>', on_vertical)
    Chainlab_2 = Label(second_frame,text="Used in the SCS models simulations and to predict proteins free energy",font="Arial 10", anchor='w', bg='gray90')
    Chainlab_2.grid(row=12, column=1, sticky="W", padx=10)
    Chainlab_2.bind('<MouseWheel>', on_vertical)
    Chainlab_3 = Label(second_frame,text="",font="Arial 10", anchor='w', bg='gray90')
    Chainlab_3.grid(row=13, column=1, sticky="W", padx=10)
    Chainlab_3.bind('<MouseWheel>', on_vertical)
    ChainEntry = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue",selectbackground="darkblue", selectforeground="white", relief="flat",highlightbackground="darkblue", highlightthickness=1)
    ChainEntry.grid(row=13, column=1, sticky='W', padx=150,pady=5)
    ChainEntry.bind('<MouseWheel>', on_vertical)

    lab10 = Label(second_frame, text="* Show running information on the screen?:", font="Arial 14", anchor='w',bg='gray90')
    lab10.grid(row=14, column=1, sticky="W", padx=5)
    lab10.bind('<MouseWheel>', on_vertical)
    lab10_2 = Label(second_frame, text="No (by default), Yes (this option increases the running time)", font="Arial 10",anchor='w', bg='gray90')
    lab10_2.grid(row=15, column=1, sticky="W", padx=10)
    lab10_2.bind('<MouseWheel>', on_vertical)
    Show_RI_menu = ttk.Combobox(second_frame, value=S_Running, state="readonly")
    Show_RI_menu.current(0)
    Show_RI_menu.grid(row=16, column=1, sticky='W', padx=150,pady=5)
    Show_RI_menu.bind('<MouseWheel>', on_vertical)

    lab7 = Label(second_frame, text=" ", font=('Arial', 14, 'bold'),bg='gray90')
    lab7.grid(row=17, column=1, sticky="W", padx=5)
    lab7.bind('<MouseWheel>', on_vertical)
    lab7_2 = Label(second_frame, text=" ## SETTINGS FOR THE SIMULATION PHASE ## ", font=('Arial', 16, 'bold'),bg='gray90')
    lab7_2.grid(row=18, column=1, sticky="W", padx=5)
    lab7_2.bind('<MouseWheel>', on_vertical)
    sep1 = Label(second_frame, text="  -------------------------------------------------------------------------------------------------------------------------------------------",font="Arial 10", anchor='w', fg='gray70', bg='gray90')
    sep1.grid(row=19, column=1, sticky="W")
    sep1.bind('<MouseWheel>', on_vertical)
    lab7_3 = Label(second_frame, text="    GENERAL SIMULATION SETTINGS ", font=('Arial', 14, 'bold'),bg='gray90')
    lab7_3.grid(row=20, column=1, sticky="W", padx=5)
    lab7_3.bind('<MouseWheel>', on_vertical)

    lab8 = Label(second_frame, text="* Number of simulations (integrer):", font="Arial 14", anchor='w', bg='gray90')
    lab8.grid(row=21, column=1, sticky="W", padx=5)
    lab8.bind('<MouseWheel>', on_vertical)
    lab8_2 = Label(second_frame, text="", font="Arial 14", anchor='w', bg='gray90')
    lab8_2.grid(row=22, column=1, sticky="W", padx=5)
    lab8_2.bind('<MouseWheel>', on_vertical)
    Numb_simu = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    Numb_simu.grid(row=22, column=1, sticky='W', padx=150, pady=5)
    Numb_simu.bind('<MouseWheel>', on_vertical)

    Oplab1 = Label(second_frame, text="* Number of processors", font="Arial 14", anchor='w', bg='gray90')
    Oplab1.grid(row=23, column=1, sticky="W", padx=5)
    Oplab1.bind('<MouseWheel>', on_vertical)
    Oplab2 = Label(second_frame, text="Max number of proccesors selected by default", font="Arial 10", anchor='w', bg='gray90')
    Oplab2.grid(row=24, column=1, sticky="W", padx=10)
    Oplab2.bind('<MouseWheel>', on_vertical)
    Oplab3 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    Oplab3.grid(row=25, column=1, sticky="W", padx=10)
    Oplab3.bind('<MouseWheel>', on_vertical)
    Numb_proce = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    Numb_proce.insert(END, os.cpu_count())
    Numb_proce.grid(row=25, column=1, sticky="W", padx=150, pady=5)
    Numb_proce.bind('<MouseWheel>', on_vertical)

    lab9 = Label(second_frame, text="* Save simulations?:", font="Arial 14", anchor='w', bg='gray90')
    lab9.grid(row=26, column=1, sticky="W", padx=5)
    lab9.bind('<MouseWheel>', on_vertical)
    lab9_2 = Label(second_frame, text="No (by default), Yes (this option requires a lot of space in the disk)", font="Arial 10",anchor='w', bg='gray90')
    lab9_2.grid(row=27, column=1, sticky="W", padx=10)
    lab9_2.bind('<MouseWheel>', on_vertical)
    Save_S_menu = ttk.Combobox(second_frame, value=S_Simulations, state="readonly")
    Save_S_menu.current(0)
    Save_S_menu.grid(row=28, column=1, sticky='W', padx=150,pady=5)
    Save_S_menu.bind('<MouseWheel>', on_vertical)

    #lab10 = Label(second_frame, text="* Show running information on the screen?:", font="Arial 14", anchor='w',bg='gray90')
    #lab10.grid(row=26, column=1, sticky="W", padx=5)
    #lab10.bind('<MouseWheel>', on_vertical)
    #lab10_2 = Label(second_frame, text="No (by default), Yes (this option increases the running time)", font="Arial 10",anchor='w', bg='gray90')
    #lab10_2.grid(row=27, column=1, sticky="W", padx=10)
    #lab10_2.bind('<MouseWheel>', on_vertical)
    #Show_RI_menu = ttk.Combobox(second_frame, value=S_Running, state="readonly")
    #Show_RI_menu.current(0)
    #Show_RI_menu.grid(row=28, column=1, sticky='W', padx=150,pady=5)
    #Show_RI_menu.bind('<MouseWheel>', on_vertical)

    sep2 = Label(second_frame, text="  -------------------------------------------------------------------------------------------------------------------------------------------",font="Arial 10", anchor='w', fg='gray70', bg='gray90')
    sep2.grid(row=29, column=1, sticky="W")
    sep2.bind('<MouseWheel>', on_vertical)
    lab11 = Label(second_frame, text="    EVOLUTIONARY HISTORY ", font=('Arial', 14, 'bold'),bg='gray90')
    lab11.grid(row=30, column=1, sticky="W", padx=5)
    lab11.bind('<MouseWheel>', on_vertical)
    lab11_2 = Label(second_frame, text="## The user should select a coalescent simulation (with user-specified parameters) or a rooted \nphylogenetic tree (provided by the user) upon which protein evolution is simulated ##", font=('Arial', 12, 'italic'),bg='gray90')
    lab11_2.grid(row=31, column=1, sticky="W", padx=5)
    lab11_2.bind('<MouseWheel>', on_vertical)

    lab12 = Label(second_frame, text="* Perform coalescent or upload phylogenetic tree?:", font="Arial 14", anchor='w',bg='gray90')
    lab12.grid(row=32, column=1, sticky="W", padx=5)
    lab12.bind('<MouseWheel>', on_vertical)

    v2 = IntVar()
    Coales = Radiobutton(second_frame, text="Coalescent", variable=v2, bg="gray90", value=0, command=lambda: Hide_General_Settings())
    Coales.bind('<MouseWheel>', on_vertical)
    Phy_tree = Radiobutton(second_frame, text="Phylogenetic tree", variable=v2, bg="gray90", value=1, command=lambda: Hide_General_Settings())
    Phy_tree.bind('<MouseWheel>', on_vertical)
    Coales.grid(column=1, row=33, sticky="W", padx=5)
    Phy_tree.grid(column=1, row=33, sticky="W", padx=150)
    lab12_2 =Label(second_frame, text="  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -", font="Arial 14", anchor='w',bg='gray90')
    lab12_2.grid(row=34, column=1, sticky="W", padx=5)
    lab12_2.bind('<MouseWheel>', on_vertical)
    
    Treename = Label(second_frame, text=" ", bg="gray90")
    lab13 = Label(second_frame, text="* Haploid or Diploid simulated data", font=('Arial', 14), anchor='w', bg='gray90')
    lab13.grid(row=35, column=1, sticky="W", padx=5)
    lab13.bind('<MouseWheel>', on_vertical)
    lab13_2 = Label(second_frame, text="", font=('Arial', 14), anchor='w', bg='gray90')
    lab13_2.grid(row=36, column=1, sticky="W", padx=5)
    lab13_2.bind('<MouseWheel>', on_vertical)
    haploid_menu = ttk.Combobox(second_frame, value=haploid, state="readonly")
    haploid_menu.current(1)
    haploid_menu.grid(row=36, column=1, sticky='W', padx=150,pady=5)
    haploid_menu.bind('<MouseWheel>', on_vertical)

    lab14 = Label(second_frame, text="* Amino acid substitution rate per site", font=('Arial', 14), anchor='w', bg='gray90')
    lab14.grid(row=37, column=1, sticky="W", padx=5)
    lab14.bind('<MouseWheel>', on_vertical)
    lab14_2 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab14_2.grid(row=38, column=1, sticky="W", padx=10)
    lab14_2.bind('<MouseWheel>', on_vertical)
    lab14_7 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab14_7.grid(row=39, column=1, sticky="W", padx=9)
    lab14_7.bind('<MouseWheel>', on_vertical)
    SR_menu = ttk.Combobox(second_frame, value=distribution, state="readonly")
    SR_menu.config(width=7)
    SR_menu.current(1)
    SR_menu.bind("<<ComboboxSelected>>", pick_SR)
    SR_menu.grid(row=39, column=1, sticky="W", padx=40,pady=5)
    SR_menu.bind('<MouseWheel>', on_vertical)
    lab14_3= Label(second_frame, text="          Example: norm              1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    lab14_3.grid(row=38, column=1, sticky="W", padx=10)
    lab14_3.bind('<MouseWheel>', on_vertical)
    lab14_4= Label(second_frame, text="", font=('Arial', 10), anchor='w', bg='gray90')
    lab14_4.grid(row=38, column=1, sticky="W", padx=10)
    lab14_4.bind('<MouseWheel>', on_vertical)
    SR_menu_1 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    SR_menu_1.config(width=12)
    SR_menu_1.grid(row=39, column=1, sticky="W", padx=140, pady=5)
    SR_menu_1.bind('<MouseWheel>', on_vertical)
    lab14_5 = Label(second_frame, text="   (t)", font=('Arial', 10), anchor='w', bg='gray90')
    lab14_5.grid(row=38, column=1, sticky="W", padx=285)
    lab14_5.bind('<MouseWheel>', on_vertical)
    SR_menu_2 = ttk.Combobox(second_frame, value=truncated, state="readonly")
    SR_menu_2.config(width=3)
    SR_menu_2.config(state="disabled")
    SR_menu_2.grid(row=39, column=1, sticky="W", padx=285, pady=5)
    SR_menu_2.bind('<MouseWheel>', on_vertical)
    lab14_6 = Label(second_frame, text="1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    lab14_6.grid(row=38, column=1, sticky="W", padx=360)
    lab14_6.bind('<MouseWheel>', on_vertical)
    SR_menu_3 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    SR_menu_3.config(width=12)
    SR_menu_3.config(state="disabled")
    SR_menu_3.grid(row=39, column=1, sticky="W", padx=360, pady=5)
    SR_menu_3.bind('<MouseWheel>', on_vertical)
    
    lab15 = Label(second_frame, text="* Effective population size (N)(integrer)", font=('Arial', 14), anchor='w', bg='gray90')
    lab15.grid(row=40, column=1, sticky="W", padx=5)
    lab15.bind('<MouseWheel>', on_vertical)
    lab15_2 = Label(second_frame, text="", font=('Arial', 14), anchor='w', bg='gray90')
    lab15_2.grid(row=41, column=1, sticky="W", padx=5)
    lab15_2.bind('<MouseWheel>', on_vertical)
    Ne = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    Ne.grid(row=41, column=1, sticky='W', padx=150,pady=5)
    Ne.bind('<MouseWheel>', on_vertical)

    Additional_parameter = ttk.Button(second_frame, text="+ Optional parameters", command=Hide_col_Settings)
    Additional_parameter.grid(row=42, column=1, sticky='W', padx=5)
    Additional_parameter.bind('<MouseWheel>', on_vertical)
    s.configure('TButton', font=('Arial', 12, 'bold'), background='gray90', relief="flat")
    s.map('TButton', background=[('active', 'gray90')])

    lab16 = Label(second_frame, text="Sampling at different times", font=('Arial', 14), anchor='w', bg='gray90')
    lab16.bind('<MouseWheel>', on_vertical)
    lab16_2 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab16_2.bind('<MouseWheel>', on_vertical)
    sampling = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue",selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    sampling.bind('<MouseWheel>', on_vertical)

    lab17 = Label(second_frame, text="Generation time", font=('Arial', 14), anchor='w', bg='gray90')
    lab17.bind('<MouseWheel>', on_vertical)
    lab17_2 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab17_2.bind('<MouseWheel>', on_vertical)
    G_time = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue",selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    G_time.bind('<MouseWheel>', on_vertical)

    lab18 = Label(second_frame, text="Exponential growth rate or Demographic periods", font=('Arial', 14), anchor='w',bg='gray90')
    lab18.bind('<MouseWheel>', on_vertical)
    lab18_2 = Label(second_frame, text="Exponential growth rate (0) or demographic periods (1)", font=('Arial', 10),anchor='w', bg='gray90')
    lab18_2.bind('<MouseWheel>', on_vertical)
    GR_menu = ttk.Combobox(second_frame, value=GrowthRate, state="readonly")
    GR_menu.current(0)
    GR_menu.bind("<<ComboboxSelected>>", pick_GR)
    GR_menu.bind('<MouseWheel>', on_vertical)
    lab19 = Label(second_frame, text="Exponential growth per individual per generation", font=('Arial', 10),anchor='w', bg='gray90')
    lab19.bind('<MouseWheel>', on_vertical)
    lab19_2 = Label(second_frame, text="Example:       1e-5", font=('Arial', 10), anchor='w', bg='gray90')
    lab19_2.bind('<MouseWheel>', on_vertical)
    GR_entry_1 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue",selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    GR_entry_1.bind('<MouseWheel>', on_vertical)
    GR_entry_1.configure(state="disabled")
    lab20 = Label(second_frame, text="Number periods + N0 + Nend + + duration per period", font=('Arial', 10), anchor='w', bg='gray90')
    lab20.bind('<MouseWheel>', on_vertical)
    lab20_2 = Label(second_frame, text="3 1000 1250 1000 1300 1550 2000 1560 1000 3000", font=('Arial', 10), anchor='w',bg='gray90')
    lab20_2.bind('<MouseWheel>', on_vertical)
    GR_entry_2 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue",selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    GR_entry_2.bind('<MouseWheel>', on_vertical)
    GR_entry_2.configure(state="disabled")

    lab21 = Label(second_frame, text="Migration model and population structure", font=('Arial', 14), anchor='w',bg='gray90')
    lab21.bind('<MouseWheel>', on_vertical)
    lab21_2 = Label(second_frame, text="See Manual for more information. Example: 2 2 3 3", font=('Arial', 10),anchor='w', bg='gray90')
    lab21_2.bind('<MouseWheel>', on_vertical)
    MM_entry = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    MM_entry.bind('<MouseWheel>', on_vertical)

    lab22 = Label(second_frame, text="Migration rate", font=('Arial', 14), anchor='w', bg='gray90')
    lab22.bind('<MouseWheel>', on_vertical)
    lab22_2 = Label(second_frame, text="See Manual for more information. Example: 3 100 800 0.002 0.001 0.003", font=('Arial', 10), anchor='w', bg='gray90')
    lab22_2.bind('<MouseWheel>', on_vertical)
    MR_entry = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue",selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    MR_entry.bind('<MouseWheel>', on_vertical)

    lab23 = Label(second_frame, text="Events of convergence of demes", font=('Arial', 14), anchor='w', bg='gray90')
    lab23.bind('<MouseWheel>', on_vertical)
    lab23_2 = Label(second_frame, text="See Manual for more information. Example: 3 1 2 400 3 4 1900 5 6 2000", font=('Arial', 10), anchor='w', bg='gray90')
    lab23_2.bind('<MouseWheel>', on_vertical)
    CD_entry = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    CD_entry.bind('<MouseWheel>', on_vertical)

    lab24 = Label(second_frame, text="User-specified phylogenetic tree/s", font=('Arial', 14), anchor='w', bg='gray90')
    lab24.bind('<MouseWheel>', on_vertical)
    lab24_2 = Label(second_frame, text="", font=('Arial', 14), anchor='w', bg='gray90')
    lab24_2.bind('<MouseWheel>', on_vertical)
    upload_PT = Button(second_frame, text="Upload ↑", command=Name_Tree, bg='white', highlightbackground='gray90')

    sep3 = Label(second_frame, text="  -------------------------------------------------------------------------------------------------------------------------------------------",font="Arial 10", anchor='w', fg='gray70', bg='gray90')
    sep3.grid(row=70, column=1, sticky="W")
    sep3.bind('<MouseWheel>', on_vertical)
    lab25 = Label(second_frame, text="    SUBSTITUTION MODEL ", font=('Arial', 14, 'bold'), bg='gray90')
    lab25.grid(row=71, column=1, sticky="W", padx=5)
    lab25.bind('<MouseWheel>', on_vertical)

    lab26 = Label(second_frame, text="* Model of amino acid substitution and template selection", font=('Arial', 14), anchor='w', bg='gray90')
    lab26.grid(row=72, column=1, sticky="W", padx=5)
    lab26.bind('<MouseWheel>', on_vertical)
    # valores
    emp_valores = StringVar()
    emp_valores.set("Blosum62 CpRev Dayhoff DayhoffDCMUT HIVb HIVw JTT JonesDCMUT LG Mtart Mtrev24 RtRev VT WAG UserEAAM")
    # crear la lista
    lab26_2 = Label(second_frame, text="Empirical substitution models", bg="gray90")
    lab26_2.grid(row=73, column=1, sticky="W", padx=5)
    lab26_2.bind('<MouseWheel>', on_vertical)
    #emp_lstbox = Listbox(second_frame, listvariable=emp_valores, selectmode=MULTIPLE, width=20,height=3)
    #emp_lstbox.grid(row=74, column=1, sticky="W", padx=10)
    #emp_lstbox.bind("<<ListboxSelect>>", pick_amino)
    emp_but = Button(second_frame, text="Select", command=openEM, bg='gray90',highlightbackground='gray90')
    emp_but.configure(width=15)
    emp_but.grid(row=74, column=1, sticky='W', padx=30, pady=10)
    emp_but.bind('<MouseWheel>', on_vertical)
    vF = IntVar()
    vN = IntVar()
    Struc_Model_F = Checkbutton(second_frame, text="Fitness", variable=vF, bg="gray90", command=checkF)
    Struc_Model_N = Checkbutton(second_frame, text="Neutral", variable=vN, bg="gray90")
    Struc_Model_F.grid(row=74, column=1, sticky="W", padx=250)
    Struc_Model_N.grid(row=74, column=1, sticky="W", padx=370)
    Struc_Model_F.select()
    Struc_Model_N.select()
    lab26_3 = Label(second_frame, text="SCS models", bg="gray90")
    lab26_3.grid(row=73, column=1, sticky="W", padx=250)
    lab26_3.bind('<MouseWheel>', on_vertical)

    lab27 = Label(second_frame, text="* Amino acid frequencies", font=('Arial', 14), anchor='w', bg='gray90')
    lab27.grid(row=75, column=1, sticky="W", padx=5)
    lab27.bind('<MouseWheel>', on_vertical)
    lab27_2 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab27_2.grid(row=76, column=1, sticky="W", padx=10)
    lab27_2.bind('<MouseWheel>', on_vertical)
    lab27_3 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab27_3.grid(row=78, column=1, sticky="W", padx=10)
    lab27_3.bind('<MouseWheel>', on_vertical)
    A_frequencies_menu=ttk.Combobox(second_frame, value=amino_distribution, width=17)
    A_frequencies_menu.bind("<<ComboboxSelected>>", pick_amino)
    A_frequencies_menu.current(0)
    A_frequencies_menu.grid(row=77, column=1, sticky='W', padx=10, pady=2)
    A_frequencies_menu.bind('<MouseWheel>', on_vertical)
    A_frequencies_menu.configure(state="disabled")
    A_frequencies = Entry(second_frame, width=18, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    A_frequencies.grid(row=78, column=1, sticky='W', padx=11)
    A_frequencies.insert(0, "0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05")
    A_frequencies.configure(state="disabled")
    A_frequencies.bind('<MouseWheel>', on_vertical)

    lab127 = Label(second_frame, text="* Temperature", font=('Arial', 12), anchor='w', bg='gray90')
    lab127.grid(row=75, column=1, sticky="W", padx=250)
    lab127.bind('<MouseWheel>', on_vertical)
    lab127_1 = Label(second_frame, text="", font=('Arial', 12), anchor='w', bg='gray90')
    lab127_1.grid(row=75, column=1, sticky="W", padx=460)
    lab127_1.bind('<MouseWheel>', on_vertical)
    Temp = Entry(second_frame, width=9, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    Temp.grid(row=75, column=1, sticky='W', padx=370)
    Temp.insert(0, "1.78")
    Temp.bind('<MouseWheel>', on_vertical)
    lab128 = Label(second_frame, text="* Unfolded entropy", font=('Arial', 12), anchor='w', bg='gray90')
    lab128.grid(row=76, column=1, sticky="W", padx=250)
    lab128.bind('<MouseWheel>', on_vertical)
    lab128_1 = Label(second_frame, text="", font=('Arial', 12), anchor='w', bg='gray90')
    lab128_1.grid(row=76, column=1, sticky="W", padx=460)
    lab128_1.bind('<MouseWheel>', on_vertical)
    S0 = Entry(second_frame, width=9, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    S0.grid(row=76, column=1, sticky='W', padx=370)
    S0.insert(0, "0.05")
    S0.bind('<MouseWheel>', on_vertical)
    lab129 = Label(second_frame, text="* Misfolded entropy", font=('Arial', 12), anchor='w', bg='gray90')
    lab129.grid(row=77, column=1, sticky="W", padx=250)
    lab129.bind('<MouseWheel>', on_vertical)
    lab129_1 = Label(second_frame, text="", font=('Arial', 12), anchor='w', bg='gray90')
    lab129_1.grid(row=77, column=1, sticky="W", padx=460)
    lab128_1.bind('<MouseWheel>', on_vertical)
    SC1 = Entry(second_frame, width=9, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    SC1.grid(row=77, column=1, sticky='W', padx=370)
    SC1.insert(0, "0.05")
    SC1.bind('<MouseWheel>', on_vertical)
    lab130 = Label(second_frame, text="* Entropy offset", font=('Arial', 12), anchor='w', bg='gray90')
    lab130.grid(row=78, column=1, sticky="W", padx=250)
    lab130.bind('<MouseWheel>', on_vertical)
    lab130_1 = Label(second_frame, text="", font=('Arial', 12), anchor='w', bg='gray90')
    lab130_1.grid(row=78, column=1, sticky="W", padx=460)
    lab130_1.bind('<MouseWheel>', on_vertical)
    SC0 = Entry(second_frame, width=9, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    SC0.grid(row=78, column=1, sticky='W', padx=370)
    SC0.insert(0, "0.0")
    SC0.bind('<MouseWheel>', on_vertical)
    lab131 = Label(second_frame, text="* REM cumulant", font=('Arial', 12), anchor='w', bg='gray90')
    lab131.grid(row=79, column=1, sticky="W", padx=250)
    lab131.bind('<MouseWheel>', on_vertical)
    lab131_1 = Label(second_frame, text="", font=('Arial', 12), anchor='w', bg='gray90')
    lab131_1.grid(row=79, column=1, sticky="W", padx=460)
    lab131_1.bind('<MouseWheel>', on_vertical)
    REM3 = Entry(second_frame, width=9, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    REM3.grid(row=79, column=1, sticky='W', padx=370)
    REM3.insert(0, "0")
    REM3.bind('<MouseWheel>', on_vertical)
    lab132 = Label(second_frame, text="* Fitness pop size", font=('Arial', 12), anchor='w', bg='gray90')
    lab132.grid(row=80, column=1, sticky="W", padx=250)
    lab132.bind('<MouseWheel>', on_vertical)
    lab132_1 = Label(second_frame, text="", font=('Arial', 12), anchor='w', bg='gray90')
    lab132_1.grid(row=80, column=1, sticky="W", padx=460)
    lab132_1.bind('<MouseWheel>', on_vertical)
    NPOP = Entry(second_frame, width=9, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    NPOP.grid(row=80, column=1, sticky='W', padx=370)
    NPOP.insert(0, "10")
    NPOP.bind('<MouseWheel>', on_vertical)

    Additional_parameters = ttk.Button(second_frame, text="+ Optional parameters", command=Hide_SM_Settings)
    Additional_parameters.grid(row=79, column=1, sticky="W")
    Additional_parameters.bind('<MouseWheel>', on_vertical)
    s.configure('TButton', font=('Arial', 12, 'bold'), background='gray90', relief="flat")
    s.map('TButton', background=[('active', 'gray90')])

    lab28 = Label(second_frame, text="Rate of heterogeneity across sites (+G)", font=('Arial', 14), anchor='w', bg='gray90')
    lab28.bind('<MouseWheel>', on_vertical)
    lab28_2 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab28_2.bind('<MouseWheel>', on_vertical)
    lab28_3 = Label(second_frame, text="Example:norm                         1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    lab28_3.bind('<MouseWheel>', on_vertical)
    G_menu = ttk.Combobox(second_frame, value=distribution, state="readonly")
    G_menu.config(width=10)
    G_menu.bind("<<ComboboxSelected>>", pick_G)
    G_menu.bind('<MouseWheel>', on_vertical)
    G_menu_1 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    G_menu_1.config(width=12)
    G_menu_1.bind('<MouseWheel>', on_vertical)
    lab28_4 = Label(second_frame, text="   (t)", font=('Arial', 10), anchor='w', bg='gray90')
    lab28_4.bind('<MouseWheel>', on_vertical)
    G_menu_2 = ttk.Combobox(second_frame, value=truncated, state="readonly")
    G_menu_2.config(width=3)
    G_menu_2.bind('<MouseWheel>', on_vertical)
    lab28_5 = Label(second_frame, text="1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    lab28_5.bind('<MouseWheel>', on_vertical)
    G_menu_3 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue",selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    G_menu_3.config(width=12)
    G_menu_3.bind('<MouseWheel>', on_vertical)

    lab29 = Label(second_frame, text="Proportion of invariable sites (+I)", font=('Arial', 14), anchor='w',bg='gray90')
    lab29.bind('<MouseWheel>', on_vertical)
    lab29_2 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab29_2.bind('<MouseWheel>', on_vertical)
    I_menu = ttk.Combobox(second_frame, value=distribution, state="readonly")
    I_menu.config(width=10)
    I_menu.bind("<<ComboboxSelected>>", pick_I)
    I_menu.bind('<MouseWheel>', on_vertical)
    lab29_3 = Label(second_frame, text="Example:norm                         1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    lab29_3.bind('<MouseWheel>', on_vertical)
    I_menu_1 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    I_menu_1.config(width=12)
    I_menu_1.bind('<MouseWheel>', on_vertical)
    lab29_4 = Label(second_frame, text="   (t)", font=('Arial', 10), anchor='w', bg='gray90')
    lab29_4.bind('<MouseWheel>', on_vertical)
    I_menu_2 = ttk.Combobox(second_frame, value=truncated, state="readonly")
    I_menu_2.config(width=3)
    I_menu_2.bind('<MouseWheel>', on_vertical)
    lab29_5 = Label(second_frame, text="1.0e-8 or 1.0e-8 1.0e-5", font=('Arial', 10), anchor='w', bg='gray90')
    lab29_5.bind('<MouseWheel>', on_vertical)
    I_menu_3 = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    I_menu_3.config(width=12)
    I_menu_3.bind('<MouseWheel>', on_vertical)

    subsection2 = Label(second_frame, text=" ## SETTINGS FOR THE ESTIMATION PHASE ## ", font=("Arial", "16", "bold"), anchor='w', bg='gray90')
    subsection2.grid(row=90, column=1, sticky="W", padx=5, pady=5)
    subsection2.bind('<MouseWheel>', on_vertical)
    sep4 = Label(second_frame, text="  -------------------------------------------------------------------------------------------------------------------------------------------",font="Arial 10", anchor='w', fg='gray70', bg='gray90')
    sep4.grid(row=91, column=1, sticky="W")
    sep4.bind('<MouseWheel>', on_vertical)

    lab30 = Label(second_frame, text="* ABC iterations", font=('Arial', 14), anchor='w', bg='gray90')
    lab30.grid(row=92, column=1, sticky="W", padx=5)
    lab30.bind('<MouseWheel>', on_vertical)
    lab30_1 = Label(second_frame, text="Iterations < NumberOfSimulations", font="Arial 10", anchor='w', bg='gray90')
    lab30_1.grid(row=93, column=1, sticky="W", padx=10)
    lab30_1.bind('<MouseWheel>', on_vertical)
    lab30_2 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab30_2.grid(row=94, column=1, sticky="W", padx=10)
    lab30_2.bind('<MouseWheel>', on_vertical)
    ABC_itera = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    ABC_itera.grid(row=94, column=1, sticky='W', padx=150, pady=5)
    ABC_itera.bind('<MouseWheel>', on_vertical)

    lab31 = Label(second_frame, text="* ABC tolerance", font=('Arial', 14), anchor='w', bg='gray90')
    lab31.grid(row=95, column=1, sticky="W", padx=5)
    lab31.bind('<MouseWheel>', on_vertical)
    lab31_1 = Label(second_frame, text="Proportion of acepted simulations. Example: 0.01", font="Arial 10", anchor='w', bg='gray90')
    lab31_1.grid(row=96, column=1, sticky="W", padx=10)
    lab31_1.bind('<MouseWheel>', on_vertical)
    lab31_2 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab31_2.grid(row=97, column=1, sticky="W", padx=10)
    lab31_2.bind('<MouseWheel>', on_vertical)
    ABC_tolerance = Entry(second_frame, width=21, bd=1, bg="white", highlightcolor="lightblue", selectbackground="darkblue", selectforeground="white", relief="flat", highlightbackground="darkblue", highlightthickness=1)
    ABC_tolerance.grid(row=97, column=1, sticky='W', padx=150, pady=5)
    ABC_tolerance.bind('<MouseWheel>', on_vertical)

    lab32 = Label(second_frame, text="* ABC method", font=('Arial', 14), anchor='w', bg='gray90')
    lab32.grid(row=98, column=1, sticky="W", padx=5)
    lab32.bind('<MouseWheel>', on_vertical)
    lab32_1 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab32_1.grid(row=99, column=1, sticky="W", padx=10)
    lab32_1.bind('<MouseWheel>', on_vertical)
    lab32_2 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab32_2.grid(row=100, column=1, sticky="W", padx=10)
    lab32_2.bind('<MouseWheel>', on_vertical)
    ABCMethod_menu = ttk.Combobox(second_frame, value=ABCMethod, state="readonly")
    ABCMethod_menu.current(0)
    ABCMethod_menu.grid(row=100, column=1, sticky='W', padx=150, pady=5)
    ABCMethod_menu.bind('<MouseWheel>', on_vertical)

    lab33 = Label(second_frame, text="* Summary statistics to use.", font=('Arial', 14), anchor='w', bg='gray90')
    lab33.grid(row=101, column=1, sticky="W", padx=5)
    lab33.bind('<MouseWheel>', on_vertical)
    lab33_1 = Label(second_frame, text="See Manual for more information ", font="Arial 10", anchor='w', bg='gray90')
    lab33_1.grid(row=102, column=1, sticky="W", padx=10)
    lab33_1.bind('<MouseWheel>', on_vertical)
    lab33_2 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab33_2.grid(row=103, column=1, sticky="W", padx=10)
    lab33_2.bind('<MouseWheel>', on_vertical)
    SubModels = Button(second_frame, text="Summary stadistics", command=openSS, bg='white', highlightbackground='gray90')
    SubModels.grid(row=103, column=1, sticky='W', padx=170, pady=5)
    SubModels.bind('<MouseWheel>', on_vertical)

    lab34 = Label(second_frame, text="* Multiple pages", font=('Arial', 12), anchor='w', bg='gray90')
    lab34.grid(row=104, column=1, sticky="W", padx=5)
    lab34.bind('<MouseWheel>', on_vertical)
    lab34_1 = Label(second_frame, text="PDF documents with multiple pages. No, Yes ", font="Arial 10", anchor='w', bg='gray90')
    lab34_1.grid(row=105, column=1, sticky="W", padx=10)
    lab34_1.bind('<MouseWheel>', on_vertical)
    lab34_2 = Label(second_frame, text="", font="Arial 10", anchor='w', bg='gray90')
    lab34_2.grid(row=106, column=1, sticky="W", padx=10)
    lab34_2.bind('<MouseWheel>', on_vertical)
    MultiPages_menu = ttk.Combobox(second_frame, value=MultiPages, state="readonly")
    MultiPages_menu.grid(row=106, column=1, sticky='W', padx=150, pady=5)
    MultiPages_menu.bind('<MouseWheel>', on_vertical)
    MultiPages_menu.current(1)

    lab43 = Label(second_frame,text="  -------------------------------------------------------------------------------------------------------------------------------------------",font="Arial 10", anchor='w', fg='gray70', bg='gray90')
    lab43.grid(row=107, column=1, sticky="W", pady=5)
    lab43.bind('<MouseWheel>', on_vertical)



    save_button = Button(second_frame, text="Save and Continue ➔", command=lambda:save_selected_values(),  bg='white', highlightbackground='gray90')
    save_button.config(font=('Arial', 14, 'bold'))
    save_button.grid(row=700, column=1, pady=5, sticky='W', padx=150)
    save_button.bind('<MouseWheel>', on_vertical)



#Splash Screen timer
splash_root.after(2000, main_window)





mainloop()
