"""
Main module for affecttracker_validation.

Author: Antonin Fourcade
Years: 2025
"""

# %% Import
import logging

from affecttracker_validation.configs import config, params, paths

# %% Set global vars & paths >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o
logger = logging.getLogger(__name__)


# %% Functions >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o


def main():
    """
    Run the main.

    This runs some example functions and logs a message, which demonstrates, how the `research-project` is structured.
    """
    # First, we extract a variable from the config object.
    # The variable 'service_x.api_key' comes from the private config file: `private_config.toml`.
    print(f"Access service x using my private key: {config.service_x.api_key}")

    # Run a placeholder function from the preprocessing module
    foo()

    # Print a placeholder path from the paths object, which was generated from the (public) config file: `config.toml`.
    print(f"{paths.results.GLM}/{params.weight_decay}/")

    # Finally, see how the pre-configured logger works
    logger.info(
        "My first log entry. Use pre-configured loggers! You can change logging.configs in 'code/configs/config.toml'"
    )


# %% __main__  >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o

if __name__ == "__main__":
    # Run main
    main()

# o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o >><< o END
