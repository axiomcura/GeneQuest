from datetime import datetime


def generate_unique_id():
    """Generates a unique ID based on the datetime format:
    %m%d%y-%H%M%S
    """
    unique_id = datetime.now().strftime("%m%d%y-%H%M%S")
    return unique_id
