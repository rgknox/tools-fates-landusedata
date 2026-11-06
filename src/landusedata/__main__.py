from landusedata._main import main

# Guard against import time side effects
if __name__ == '__main__':
    raise SystemExit(main())
