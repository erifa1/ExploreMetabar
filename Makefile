# Release automation for ExploreMetabar.
#
#   make release VERSION=3.1.0
#
# Bumps the version, tags the package repo (forge.inrae.fr/umrf/exploremetabar),
# then syncs + pushes the SK8 deploy repo (sk8/sk8-apps/tbi/explore-metabar),
# which fires the SK8 build -> auto-deploy. See the deploy contract in CLAUDE.md.
#
# DEPLOY_REPO is the local checkout of the SK8 deploy repo (override per machine).

DEPLOY_REPO ?= $(HOME)/git-repo/explore-metabar
VERSION ?=

.PHONY: release sync-deploy

## Bump version, tag repo A, push, then sync the deploy repo.
release:
	@test -n "$(VERSION)" || { echo "ERROR: VERSION= is required (e.g. make release VERSION=3.1.0)"; exit 1; }
	@test "$$(git rev-parse --abbrev-ref HEAD)" = "master" || { echo "ERROR: release must be cut from master"; exit 1; }
	@git diff --quiet && git diff --cached --quiet || { echo "ERROR: repo has uncommitted changes"; exit 1; }
	sed -i 's/^Version: .*/Version: $(VERSION)/' DESCRIPTION
	sed -i 's/^  golem_version: .*/  golem_version: $(VERSION)/' inst/golem-config.yml
	@grep -q "ExploreMetabar $(VERSION)" NEWS.md || echo ">> WARNING: NEWS.md has no section for $(VERSION) -- add one before releasing."
	@git diff --quiet || git commit -am "Release $(VERSION)"
	git tag -a v$(VERSION) -m "ExploreMetabar $(VERSION)"
	git push origin master --follow-tags
	$(MAKE) sync-deploy VERSION=$(VERSION)

## Copy the tagged renv.lock + pin ARG EM_REF into the deploy repo, then push (-> SK8 build).
sync-deploy:
	@test -n "$(VERSION)" || { echo "ERROR: VERSION= is required"; exit 1; }
	@test -d "$(DEPLOY_REPO)/.git" || { echo "ERROR: DEPLOY_REPO is not a git repo: $(DEPLOY_REPO)"; exit 1; }
	git -C "$(DEPLOY_REPO)" checkout main
	git -C "$(DEPLOY_REPO)" pull --ff-only
	cp renv.lock "$(DEPLOY_REPO)/renv.lock"
	sed -i 's/^ARG EM_REF=.*/ARG EM_REF=v$(VERSION)/' "$(DEPLOY_REPO)/Dockerfile"
	git -C "$(DEPLOY_REPO)" add renv.lock Dockerfile
	@git -C "$(DEPLOY_REPO)" diff --cached --quiet || git -C "$(DEPLOY_REPO)" commit -m "Deploy ExploreMetabar v$(VERSION)"
	git -C "$(DEPLOY_REPO)" push origin main
