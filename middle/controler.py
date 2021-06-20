from back.process_data import Program
import eel


program = Program()


@eel.expose
def select_file(type):
    return program.get_data(type)


@eel.expose
def run_programm(length=None, threshold=None):
    return program.run(length, threshold)


@eel.expose
def open_result():
    return program.open_file()
