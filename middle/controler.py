from back.process_data import Program
import eel


program = Program()


@eel.expose
def select_file(type):
    return program.get_data(type)


@eel.expose
def run_programm(threshold=None):
    return program.run(threshold)


@eel.expose
def open_result():
    return program.open_file()

@eel.expose
def open_sliding(threshold=None):
    return program.sliding(threshold)


@eel.expose
def show_plot():
    return program.show_graphic()
