import sys
import time


def get_current_time() -> str:
    return time.strftime('%H:%M:%S', time.localtime())


def write_direct_message(message: str):
    curr_time_str = get_current_time()
    sys.stdout.write(f'{curr_time_str} --- {message}\n')
    sys.stdout.flush()


def debug(message: str):
    write_direct_message(f'DEBUG: {message}')


def info(message: str):
    write_direct_message(f'INFO: {message}')


def write_direct_message_err(message: str):
    curr_time_str = get_current_time()
    sys.stderr.write(f'{curr_time_str} --- {message}\n')
    sys.stderr.flush()


def warning(message: str):
    write_direct_message_err(f'WARNING: {message}')


def error(message: str):
    write_direct_message_err(f'ERROR: {message}')


def critical(message: str):
    write_direct_message_err(f'CRITICAL: {message}')


__all__ = ['debug', 'info', 'warning', 'error', 'critical']
