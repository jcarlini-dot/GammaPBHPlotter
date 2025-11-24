def test_import_and_version():
    import gammapbh as gp
    assert hasattr(gp, "__version__")
    assert isinstance(gp.__version__, str)