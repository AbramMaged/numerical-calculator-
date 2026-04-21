# gui.py
import customtkinter as ctk
import pandas as pd
from fpdf import FPDF
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter.messagebox as messagebox
import sympy as sp
# Import logic from main.py
import main as backend

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("dark-blue")

# Custom Colors
BG_COLOR = "#0f071f"
PANEL_COLOR = "#1b0f2e"
DARK_PANEL_COLOR = "#160a26"
PINK_ACCENT = "#d1208a"
PURPLE_ACCENT = "#8a2be2"
TEXT_COLOR = "#e5d0f5"

class NumericalApp(ctk.CTk): #
    def __init__(self):
        super().__init__()
        # main window
        self.title("Numerical Methods Solver")
        self.geometry("1100x750")
        self.configure(fg_color=BG_COLOR)
        
        self.iterations_data = [] 
        self.final_root = "-"
        self.final_error = "-"
        self.final_status = "-"
        self.eval_func = None
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.fig.patch.set_facecolor(DARK_PANEL_COLOR)
        self.ax.set_facecolor(DARK_PANEL_COLOR)
        self.ax.tick_params(colors=TEXT_COLOR)
        for spine in self.ax.spines.values():
            spine.set_color(PURPLE_ACCENT)
        self.canvas = None
        
        self._build_ui()
        self.update_input_fields_by_method("Bisection")
        self.bind('<Return>', self.solve_root)

    # ui header
    def _build_ui(self):
        header_frame = ctk.CTkFrame(self, fg_color=BG_COLOR, corner_radius=0)
        header_frame.pack(fill="x", pady=10, padx=15)
        #tle3 3yne 3lshan 2gyb 2l symbol da ⚙
        title_lbl = ctk.CTkLabel(header_frame, text="⚙📠 Numerical Methods Solver", font=ctk.CTkFont(size=20, weight="bold"), text_color=PINK_ACCENT)
        title_lbl.pack(side="left")
        
        self.tabview = ctk.CTkTabview(self, fg_color="transparent", text_color=TEXT_COLOR, segmented_button_selected_color=PURPLE_ACCENT, segmented_button_selected_hover_color=PINK_ACCENT)
        self.tabview.pack(fill="both", expand=True, padx=15, pady=5)
        self.tabview.add("Root Finding")
        self.tabview.add("Linear Systems")
        
        self._build_root_finding_tab()
        self._build_linear_systems_tab()
    #SIDE TAP    
    def _build_root_finding_tab(self):
        main_frame = ctk.CTkFrame(self.tabview.tab("Root Finding"), fg_color="transparent")
        main_frame.pack(fill="both", expand=True)

        left_col = ctk.CTkFrame(main_frame, fg_color="transparent", width=400)
        left_col.pack(side="left", fill="y", expand=False)
        
        right_col = ctk.CTkFrame(main_frame, fg_color=DARK_PANEL_COLOR, corner_radius=10)
        right_col.pack(side="right", fill="both", expand=True, padx=(15,0))
        
        method_lbl = ctk.CTkLabel(left_col, text="Method", font=ctk.CTkFont(size=12), text_color=TEXT_COLOR)
        method_lbl.pack(anchor="w", pady=(0, 2))
        self.method_var = ctk.StringVar(value="Bisection")
        method_menu = ctk.CTkOptionMenu(left_col, variable=self.method_var, values=["Bisection", "False Position", "Fixed Point", "Newton-Raphson", "Secant"], fg_color=PURPLE_ACCENT, button_color=PINK_ACCENT, button_hover_color=PURPLE_ACCENT, command=self.update_input_fields_by_method)
        method_menu.pack(fill="x", pady=(0, 15))
        
        inputs_panel = ctk.CTkFrame(left_col, fg_color=PANEL_COLOR, corner_radius=10)
        inputs_panel.pack(fill="x", pady=5)
        lbl = ctk.CTkLabel(inputs_panel, text="⚙ Input Parameters", font=ctk.CTkFont(size=14, weight="bold"), text_color=PINK_ACCENT)
        lbl.pack(anchor="w", padx=10, pady=10)
        
        ctk.CTkLabel(inputs_panel, text="f(x) =", text_color=TEXT_COLOR).pack(anchor="w", padx=10)
        #e.g to help the user with the equation input-----------------------------------------------------------------#
        self.entry_fx = ctk.CTkEntry(inputs_panel, placeholder_text="e.g. x**3 - x - 2", fg_color=BG_COLOR, border_color=PURPLE_ACCENT)
        self.entry_fx.pack(fill="x", padx=10, pady=(0, 10))
        
        #xl and xu ui
        ab_frame = ctk.CTkFrame(inputs_panel, fg_color="transparent")
        ab_frame.pack(fill="x", padx=10, pady=(0, 10))
        self.lbl_a = ctk.CTkLabel(ab_frame, text="xl =", text_color=TEXT_COLOR)
        self.lbl_a.grid(row=0, column=0, sticky="w", padx=(0,5))
        self.entry_a = ctk.CTkEntry(ab_frame, width=80, fg_color=BG_COLOR, border_color=PURPLE_ACCENT)
        self.entry_a.grid(row=1, column=0, sticky="w", padx=(0,15))
        self.entry_a.insert(0, "1") #e.g xl-----------------------------------------------------------------#

        self.lbl_b = ctk.CTkLabel(ab_frame, text="xu =", text_color=TEXT_COLOR)
        self.lbl_b.grid(row=0, column=1, sticky="w", padx=(0,5))
        self.entry_b = ctk.CTkEntry(ab_frame, width=80, fg_color=BG_COLOR, border_color=PURPLE_ACCENT)
        self.entry_b.grid(row=1, column=1, sticky="w")
        self.entry_b.insert(0, "2") #e.g xu-----------------------------------------------------------------#
        
        tol_frame = ctk.CTkFrame(inputs_panel, fg_color="transparent")
        tol_frame.pack(fill="x", padx=10, pady=(0, 15))
        # tolarance ui
        ctk.CTkLabel(tol_frame, text="Tolerance", text_color=TEXT_COLOR).grid(row=0, column=0, sticky="w")
        self.entry_tol = ctk.CTkEntry(tol_frame, width=80, fg_color=BG_COLOR, border_color=PURPLE_ACCENT)
        self.entry_tol.grid(row=1, column=0, sticky="w", padx=(0,10))
        self.entry_tol.insert(0, ".5") #e.g tolerance-----------------------------------------------------------------#
        
        #max iterations ui
        ctk.CTkLabel(tol_frame, text="Max Iterations", text_color=TEXT_COLOR).grid(row=0, column=1, sticky="w")
        self.entry_maxit = ctk.CTkEntry(tol_frame, width=80, fg_color=BG_COLOR, border_color=PURPLE_ACCENT)
        self.entry_maxit.grid(row=1, column=1, sticky="w", padx=(0,10))
        self.entry_maxit.insert(0, "50") #e.g max iterations-----------------------------------------------------------------#
        
        #precision ui
        ctk.CTkLabel(tol_frame, text="Precision", text_color=TEXT_COLOR).grid(row=0, column=2, sticky="w")
        self.entry_prec = ctk.CTkEntry(tol_frame, width=60, fg_color=BG_COLOR, border_color=PURPLE_ACCENT)
        self.entry_prec.grid(row=1, column=2, sticky="w")
        self.entry_prec.insert(0, "4") #e.g precision-----------------------------------------------------------------#
        
        #buttons ui
        btns_frame = ctk.CTkFrame(left_col, fg_color="transparent")
        btns_frame.pack(fill="x", pady=10)
        ctk.CTkButton(btns_frame, text="Solve", fg_color=PINK_ACCENT, hover_color="#b81c7a", width=80, command=self.solve_root).pack(side="left", padx=(0,5))
        ctk.CTkButton(btns_frame, text="PDF", fg_color="transparent", border_width=1, border_color=PURPLE_ACCENT, text_color=TEXT_COLOR, width=50, command=self.export_pdf).pack(side="left", padx=5)
        ctk.CTkButton(btns_frame, text="Excel", fg_color="transparent", border_width=1, border_color=PURPLE_ACCENT, text_color=TEXT_COLOR, width=50, command=self.export_excel).pack(side="left", padx=5)
        
        #results ui
        res_panel = ctk.CTkFrame(left_col, fg_color=PANEL_COLOR, corner_radius=10)
        res_panel.pack(fill="x", pady=10)
        ctk.CTkLabel(res_panel, text="📋 Results", font=ctk.CTkFont(size=14, weight="bold"), text_color=PINK_ACCENT).pack(anchor="w", padx=10, pady=(10, 5))
        res_grid = ctk.CTkFrame(res_panel, fg_color="transparent")
        res_grid.pack(fill="x", padx=10, pady=(0, 10))
        res_grid.columnconfigure((0,1), weight=1)
        
        #results ui
        def make_res_card(parent, title, r, c):
            card = ctk.CTkFrame(parent, fg_color=DARK_PANEL_COLOR, corner_radius=6)
            card.grid(row=r, column=c, padx=5, pady=5, sticky="nsew")
            ctk.CTkLabel(card, text=title, text_color=PURPLE_ACCENT, font=ctk.CTkFont(weight="bold")).pack(anchor="w", padx=5, pady=(5,0))
            val_lbl = ctk.CTkLabel(card, text="-", text_color="white")
            val_lbl.pack(anchor="w", padx=5, pady=(0,5))
            return val_lbl
        
        self.lbl_res_root = make_res_card(res_grid, "Root", 0, 0)
        self.lbl_res_iter = make_res_card(res_grid, "Iterations", 0, 1)
        self.lbl_res_err = make_res_card(res_grid, "Error", 1, 0)
        self.lbl_res_status = make_res_card(res_grid, "Status", 1, 1)
        
        ctk.CTkLabel(right_col, text="📈 Visualization", font=ctk.CTkFont(size=14, weight="bold"), text_color=PINK_ACCENT).pack(anchor="w", padx=15, pady=10)
        self.canvas_frame = ctk.CTkFrame(right_col, fg_color="transparent")
        self.canvas_frame.pack(fill="both", expand=True, padx=10, pady=(0,10))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
    
    #linear systems tab ui
    def _build_linear_systems_tab(self):
        lin_frame = ctk.CTkFrame(self.tabview.tab("Linear Systems"), fg_color="transparent")
        lin_frame.pack(fill="both", expand=True)
        
        #control panel ui
        control_panel = ctk.CTkFrame(lin_frame, fg_color=PANEL_COLOR, corner_radius=10)
        control_panel.pack(side="left", fill="y", padx=10, pady=10)
        
        #method ui
        ctk.CTkLabel(control_panel, text="Method", font=ctk.CTkFont(size=12), text_color=TEXT_COLOR).pack(anchor="w", padx=15, pady=(15, 2))
        self.lin_method_var = ctk.StringVar(value="Gauss-Jordan Elimination")
        lin_menu = ctk.CTkOptionMenu(control_panel, variable=self.lin_method_var, values=["Gauss-Jordan Elimination", "LU Decomposition", "Cramer's Rule"], fg_color=PURPLE_ACCENT, button_color=PINK_ACCENT, button_hover_color=PURPLE_ACCENT)
        lin_menu.pack(fill="x", padx=15, pady=(0, 15))
        
        #matrix ui
        ctk.CTkLabel(control_panel, text="Enter Augmented Matrix (3x4)", font=ctk.CTkFont(weight="bold"), text_color=PINK_ACCENT).pack(pady=10)
        
        self.matrix_entries = []
        grid_frame = ctk.CTkFrame(control_panel, fg_color="transparent")
        grid_frame.pack(padx=15, pady=10)
        
        for i in range(3):
            row_entries = []
            for j in range(4):
                e = ctk.CTkEntry(grid_frame, width=50, fg_color=BG_COLOR, border_color=PURPLE_ACCENT)
                e.grid(row=i, column=j, padx=5, pady=5)
                e.insert(0, "0")
                if j == 3: e.configure(border_color=PINK_ACCENT)
                row_entries.append(e)
            self.matrix_entries.append(row_entries)
            
        ctk.CTkButton(control_panel, text="Solve System", fg_color=PINK_ACCENT, hover_color="#b81c7a", command=self.solve_linear).pack(pady=20, padx=15, fill="x")
        
        lin_res_panel = ctk.CTkFrame(lin_frame, fg_color=DARK_PANEL_COLOR, corner_radius=10)
        lin_res_panel.pack(side="right", fill="both", expand=True, padx=10, pady=10)
        
        ctk.CTkLabel(lin_res_panel, text="Solution (X Vector)", font=ctk.CTkFont(size=18, weight="bold"), text_color=PINK_ACCENT).pack(pady=20)
        
        self.lbl_x1 = ctk.CTkLabel(lin_res_panel, text="X1 = -", font=ctk.CTkFont(size=14), text_color=TEXT_COLOR)
        self.lbl_x1.pack(pady=5)
        self.lbl_x2 = ctk.CTkLabel(lin_res_panel, text="X2 = -", font=ctk.CTkFont(size=14), text_color=TEXT_COLOR)
        self.lbl_x2.pack(pady=5)
        self.lbl_x3 = ctk.CTkLabel(lin_res_panel, text="X3 = -", font=ctk.CTkFont(size=14), text_color=TEXT_COLOR)
        self.lbl_x3.pack(pady=5)

    def update_input_fields_by_method(self, method_name):
        if method_name in ["Newton-Raphson", "Fixed Point"]:
            self.lbl_a.configure(text="x0 =")
            self.lbl_b.grid_remove()
            self.entry_b.grid_remove() # hide them till needed
        elif method_name == "Secant":
            self.lbl_a.configure(text="x0 =")
            self.lbl_b.configure(text="x1 =")
            self.lbl_b.grid()
            self.entry_b.grid()
        else: 
            self.lbl_a.configure(text="xl =")
            self.lbl_b.configure(text="xu =")
            self.lbl_b.grid()
            self.entry_b.grid()
    
    def _parse_equation(self, eq_str):
        try:
            eq_str = eq_str.replace('^', '**')
            x = sp.Symbol('x')
            expr = sp.sympify(eq_str)
            f = sp.lambdify(x, expr, 'numpy')
            f_dash_expr = sp.diff(expr, x)
            f_dash = sp.lambdify(x, f_dash_expr, 'numpy')
            return f, f_dash
        except Exception as e:
            messagebox.showerror("Parse Error", f"Could not parse the equation: {e}")
            return None, None
    
    def solve_root(self, event=None):
        eq_str = self.entry_fx.get()
        if not eq_str:
            messagebox.showerror("Error", "Please enter an equation.")
            return
            
        f, f_dash = self._parse_equation(eq_str)
        if f is None: return
        self.eval_func = f

        method = self.method_var.get()
        try:
            tol = float(self.entry_tol.get())
            max_iter = int(self.entry_maxit.get())
            prec = int(self.entry_prec.get())
        except ValueError:
            messagebox.showerror("Input Error", "Tolerance, Max Iterations, and Precision must be numbers.")
            return

        # Interface with backend from main.py
        try:
            if method == "Bisection":
                xl, xu = float(self.entry_a.get()), float(self.entry_b.get())
                root, error, iterations, status, data = backend.run_bisection(f, xl, xu, tol, max_iter)
            elif method == "False Position":
                xl, xu = float(self.entry_a.get()), float(self.entry_b.get())
                root, error, iterations, status, data = backend.run_false_position(f, xl, xu, tol, max_iter)
            elif method == "Fixed Point":
                x0 = float(self.entry_a.get())
                root, error, iterations, status, data = backend.run_fixed_point(f, x0, tol, max_iter)
            elif method == "Newton-Raphson":
                x0 = float(self.entry_a.get())
                root, error, iterations, status, data = backend.run_newton(f, f_dash, x0, tol, max_iter)
            elif method == "Secant":
                x0, x1 = float(self.entry_a.get()), float(self.entry_b.get())
                root, error, iterations, status, data = backend.run_secant(f, x0, x1, tol, max_iter)
                
            self.iterations_data = data
            
            # =====================================================================#
            # root/error/iterations/status
            # =====================================================================#

            # 3. UPDATE RESULTS
            if root is not None:
                self.final_root = f"{root:.{prec}f}"
                if error > 0 and error < 10**(-prec):
                    self.final_error = f"{error:.{prec}e}%"
                else:
                    self.final_error = f"{error:.{prec}f}%"
                self.final_status = status
                
                self.lbl_res_root.configure(text=self.final_root)
                self.lbl_res_iter.configure(text=str(iterations))
                self.lbl_res_err.configure(text=self.final_error)
                self.lbl_res_status.configure(text=status)
                
                self._update_plot(root)
            else:
                self.lbl_res_status.configure(text="Init Bounds Error")
                
        except Exception as e:
            messagebox.showerror("Execution Error", str(e))

    # 4. UPDATE PLOT
    def _update_plot(self, root):
        self.ax.clear()
        self.ax.grid(color=PANEL_COLOR, linestyle='--', linewidth=0.5)
        self.ax.axhline(0, color=TEXT_COLOR, linewidth=1)
        self.ax.axvline(0, color=TEXT_COLOR, linewidth=1)
        if self.eval_func:
            x_vals = np.linspace(root - 5, root + 5, 400)
            y_vals = []
            for xv in x_vals:
                try: y_vals.append(self.eval_func(xv))
                except: y_vals.append(np.nan)
            self.ax.plot(x_vals, y_vals, color=PINK_ACCENT, label="f(x)")
            self.ax.plot(root, 0, 'ro', label=f"Root: {root:.4f}", markerfacecolor='yellow')
            self.ax.legend(loc="upper right", frameon=False, labelcolor=TEXT_COLOR)
        self.canvas.draw()

    def solve_linear(self):
        a = []
        try:
            for i in range(3):
                row = []
                for j in range(4):
                    row.append(float(self.matrix_entries[i][j].get()))
                a.append(row)
        except ValueError:
            messagebox.showerror("Input Error", "All matrix elements must be numbers.")
            return

        method = self.lin_method_var.get()
        x = None
        # Interface with backend from main.py
        try:
            if method == "Gauss-Jordan Elimination":
                x = backend.gauss_jordan_elimination(a)
            elif method == "LU Decomposition":
                x = backend.lu_decomposition(a)
            elif method == "Cramer's Rule":
                x = backend.cramers_rule(a)
                
            if x is not None:
                self.lbl_x1.configure(text=f"X1 = {x[0]:.6f}")
                self.lbl_x2.configure(text=f"X2 = {x[1]:.6f}")
                self.lbl_x3.configure(text=f"X3 = {x[2]:.6f}")
        except Exception as e:
            messagebox.showerror("Execution Error", str(e))

    def export_excel(self):
        if not self.iterations_data: return
        pd.DataFrame(self.iterations_data).to_excel("results_export.xlsx", index=False)
        messagebox.showinfo("Export Success", "Exported successfully to results_export.xlsx")
# ==================================================================#
# EXTRA GRADES PDF EXPORTING FUNCTION
#===================================================================#
    def export_pdf(self):
        if not self.iterations_data: return
        pdf = FPDF()
        pdf.add_page()
        pdf.set_font("Arial", size=12)
        pdf.cell(200, 10, txt="Numerical Methods Solver", ln=1, align='C')
        pdf.cell(200, 10, txt=f"Method: {self.method_var.get()}", ln=1, align='L')
        pdf.cell(200, 10, txt=f"Final Root: {self.final_root}  | Error: {self.final_error}", ln=1)
        df = pd.DataFrame(self.iterations_data)
        pdf.set_font("Courier", size=8)
        pdf.ln(10)
        pdf.cell(0, 10, txt=" | ".join(df.columns), ln=1)
        for _, row in df.iterrows():
            pdf.cell(0, 5, txt=" | ".join([f"{str(v)[:8]}" if isinstance(v, float) else str(v) for v in row]), ln=1)
        pdf.output("results_export.pdf")
        messagebox.showinfo("Export", "Exported successfully to results_export.pdf")

if __name__ == "__main__":
    app = NumericalApp()
    app.mainloop()
