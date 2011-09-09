from distutils.core import setup

PKG_NAME = "repixel"
PKG_VERSION = "1.4"
PKG_DESC = "python/numarray repixeling package"
PKG_AUTHOR = "Niklas Karlsson"
PKG_AUTHMAIL = "niklas@slac.stanford.edu"
PKG_URL = "http://www.slac.stanford.edu/~niklas/Public/"

setup(name = PKG_NAME,
      version = PKG_VERSION,
      description = PKG_DESC,
      author = PKG_AUTHOR,
      author_email = PKG_AUTHMAIL,
      url = PKG_URL,
      packages = ["repixel"],
      package_dir = {"repixel": "lib"})
