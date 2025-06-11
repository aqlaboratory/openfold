#!/usr/bin/env python3

def test_openmm():
    """Test OpenMM installation and functionality"""
    print("=" * 50)
    print("Testing OpenMM and PDBFixer")
    print("=" * 50)
    
    try:
        import openmm
        print(f"✅ OpenMM imported successfully")
        print(f"   Version: {openmm.__version__}")
        
        print(f"\nAvailable OpenMM platforms:")
        for i in range(openmm.Platform.getNumPlatforms()):
            platform = openmm.Platform.getPlatform(i)
            print(f"   {i}: {platform.getName()}")
            
    except Exception as e:
        print(f"❌ Failed to import/test OpenMM: {e}")
        return False
    
    try:
        import pdbfixer
        print(f"✅ PDBFixer imported successfully")
        
        # Test basic PDBFixer functionality
        from pdbfixer import PDBFixer
        print(f"✅ PDBFixer class imported successfully")
        
    except Exception as e:
        print(f"❌ Failed to import/test PDBFixer: {e}")
        return False
    
    print(f"\n✅ All tests passed!")
    return True

if __name__ == "__main__":
    test_openmm() 