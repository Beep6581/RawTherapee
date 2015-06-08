@echo off
for /f "tokens=*" %%a in ('hg branch') do @set hgBranch=%%a
for /f "tokens=*" %%a in ('hg parents --template "{latesttag}"') do @set hgLatesttag=%%a
for /f "tokens=*" %%a in ('hg parents --template "{latesttagdistance}"') do @set hgLatesttagdistance=%%a
for /f "tokens=*" %%a in ('hg parents --template "{node|short}"') do @set hgChangeset=%%a

echo set(HG_BRANCH %hgBranch%) > ReleaseInfo.cmake
echo set(HG_VERSION %hgLatesttag%.%hgLatesttagdistance%) >> ReleaseInfo.cmake
echo set(HG_CHANGESET %hgChangeset%) >> ReleaseInfo.cmake
echo set(HG_TAGDISTANCE %hgLatesttagdistance%) >> ReleaseInfo.cmake

