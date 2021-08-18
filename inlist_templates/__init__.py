from . import simple_star
from . import simple_g_mode
from . import simple_p_mode
from . import mixing_star

templates = [
    simple_star,
    mixing_star,
    simple_g_mode,
    simple_p_mode
]


def get_template(code: str, name: str):
    for template in templates:
        if template.code == code and template.name == name:
            return template
