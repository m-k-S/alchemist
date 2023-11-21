import subprocess
from sshconf import read_ssh_config
import webbrowser
from os.path import expanduser
import argparse


def main():
    parser = argparse.ArgumentParser(prog="Alchemist")
    subparsers = parser.add_subparsers(dest="command")

    # Define a subcommand
    start_parser = subparsers.add_parser("start")
    start_parser.add_argument("--port", help="port to open jupyter server on", default=8080)
    start_parser.add_argument("--token", help="jupyter token", default="alchemst")

    connect_parser = subparsers.add_parser("connect")
    connect_parser.add_argument("host", type=str, help="ssh hostname of remote server")
    connect_parser.add_argument("--port", help="port to connect to remote jupyter server", default="8080")

    # Parse arguments
    args = parser.parse_args()

    if args.command == "start":
        launch_server(args)

    if args.command == "connect":
        connect(args)

def launch_server(args):
    port = args.port
    token = args.token
    command = ["jupyter", "notebook", "--port={port}".format(port=port), "--NotebookApp.token={token}".format(token=token), "--no-browser"]
    subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def connect(args):
    port = args.port
    host = args.host
    c = read_ssh_config(expanduser("~/.ssh/config"))
    if host not in c.hosts():
        raise ValueError("Host not found in ssh config file")
    else: 
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW  

        command = ['ssh', '-L', '{port}:localhost:{port}'.format(port=port), host]
        subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, startupinfo=startupinfo)
        webbrowser.open('localhost:{port}'.format(port=port), new=2)

if __name__ == "__main__":
    main()
