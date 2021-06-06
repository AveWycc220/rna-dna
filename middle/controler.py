from back.process_data import Program
import eel


program = Program()


@eel.expose
def select_file(type):
    return program.get_data(type)


@eel.expose
def run_programm():
    return program.run()
