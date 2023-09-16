#------------------------------------------------------------
# Common code shared among several build management scripts
# for getting lists of
# -platform directory items
# -full path names of platform directory items
# -linknames for promoted builds (releases)
# -build names of promoted builds
#
# To employ in a program place
#  'import build_mgmt' at beginning of source.
#------------------------------------------------------------

import os

platforms = ['Linux_x86_64_intel']
#platforms = ['Linux_x86_64_intel', 'Linux_x86_64_gfortran', 'packages_intel', 'packages_gfortran']
promotion_labels = ['devel', 'current']
release_prefix = 'cesr'
builds_basedir = '/nfs/acc/libs/'

diritems = {}
for platform in platforms:
    diritems[platform] = []
    
fulldiritems = {}
for platform in platforms:
    fulldiritems[platform] = []
    
prolinks = {}
for platform in platforms:
    prolinks[platform] = []

builds = {}
for platform in platforms:
    builds[platform] = []

releases = {}
for platform in platforms:
    releases[platform] = []


active_releases = {}
for platform in platforms:
    active_releases[platform] = {}
    for promotion_label in promotion_labels:
        active_releases[platform][promotion_label] = ''


for platform in platforms:
    subdirs = []
    protected_releases = []
    basedir = builds_basedir + platform
    diritems[platform] = os.listdir(basedir)
    diritems[platform].sort()
    for item in diritems[platform]:
        item = basedir + '/' + item
        fulldiritems[platform].append(item)
    subdirs = os.listdir(basedir)
    subdirs.sort()
    for item in fulldiritems:
        if release_prefix in item:
            builds[platform].append(item)

    # Get all promoted release names
    for fulldiritem in fulldiritems[platform]:
        for promotion_label in promotion_labels:
            truename = ''
            if promotion_label in fulldiritem:
                prolinks[platform].append(os.path.split(fulldiritem)[1])
                truename = os.readlink(fulldiritem)
                releases[platform].append(truename)
            if os.path.split(fulldiritem)[1] == promotion_label:
                active_releases[platform][promotion_label] = truename


# These methods provide access to information about the release archive contents
# harvested by the above code or obtained by their own execution.

def Platforms():
    """Returns a list of all supported platforms."""
    return platforms

def DirItems():
    """Returns a dictionary of all library archive platform
    directory items (short names) for each supported platform."""
    return diritems

def FullDirItems():
    """Returns a dictionary of all library archive platform
    directory items with full pathnames for each supported platform."""
    return fulldiritems

def Prolinks(promotion_label):
    """Returns a dictionary of all promoted build (release)
    symlink names for each supported platform."""
    return prolinks

def NewLinkName(promotion_label):
    """Return the name of the next link name to be used when archiving
    a previously active release directory for the given promotion label.
    This picks the highest numerical portion of the promotion label
    type found in all supported platform directories and then adds 1 to
    it to compose the new archival promotion label.
    Ex.
      Highest old archival release label in platform dir 1 = 'devel_123'
      Highest old archival release label in platform dir 2 = 'devel_173'
       (This mismatch is unlikely to happen, but it is supported.)
      New link name returned = 'devel_174' and can be used for
      archival rotation in all supported platform directories.
      """
    latestlinks = {}
    highval = 0
    for platform in platforms:
        latestlinks[platform] = ''
        for link in prolinks[platform]:
            if promotion_label in link and '_' in link and \
                   link > latestlinks[platform]:
                latestlinks[platform] = link
        try:
            value = int(latestlinks[platform].split('_')[1])+1
            if value > highval:
                highval = value
        except IndexError:
            pass
    return promotion_label+'_'+str(highval)

def Builds():
    """Returns a dictionary of all archived builds, whether or
    not they have ever been promoted to a release, for each
    supported platform."""
    return builds

def Releases():
    """Returns dictionary of all releases, active and past,
    for each supported platform."""
    return releases

def ActiveRelNames(promotion_label):
    """Return dictionary of true release names for the given
    promotion label (typically 'current' or 'devel', defined
    above) for each supported platform."""
    names = {}
    for platform in platforms:
        names[platform] = active_releases[platform][promotion_label]
    return names

def RelNameConsistency(promotion_label):
    """Returns True if the true release name for the given
    promotion label is identical across all supported platforms,
    returns False otherwise."""
    names = []
    for platform in platforms:
        names.append(active_releases[platform][promotion_label])
    if len(set(names)) == 1:
        return True
    else:
        return False

def BuildExists(buildname):
    """Returns True if the specified build name is present
    in all supported platform archive directories, False
    otherwise."""
    for platform in platforms:
        if not os.path.exists(builds_basedir+'/'+platform+'/'+buildname):
            return False
    return True
                      
